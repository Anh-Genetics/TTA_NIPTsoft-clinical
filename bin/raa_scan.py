#!/usr/bin/env python3
"""
raa_scan.py
===========
EN: Rare Autosomal Aneuploidy (RAA) scan module.
VI: Module quét lệch bội NST thường hiếm gặp (RAA - Rare Autosomal Aneuploidy).

EN: Scans specified autosomes for whole-chromosome aneuploidies
    using stricter thresholds than the primary T21/T18/T13 module.
    Results are flagged with explicit "insufficient evidence" when
    confidence is low. Conservative no-call behavior is mandatory.
VI: Quét các NST thường được chỉ định tìm lệch bội toàn bộ NST
    bằng ngưỡng nghiêm ngặt hơn module chính T21/T18/T13.
    Kết quả được đánh dấu rõ ràng là "không đủ bằng chứng" khi
    độ tin cậy thấp. Hành vi thận trọng khi không gọi là bắt buộc.
"""

import argparse
import json
import sys
from datetime import datetime, timezone

try:
    import pandas as pd
    import numpy as np
    import pyarrow.parquet as pq
except ImportError as e:
    # EN: Fatal — required libraries missing
    # VI: Lỗi nghiêm trọng — thiếu thư viện bắt buộc
    print(f"ERROR / LỖI: Missing dependency / Thiếu thư viện: {e}", file=sys.stderr)
    sys.exit(1)


# EN: Base normalizer autosomes (T21/T18/T13 excluded to avoid confounding)
# VI: Danh sách NST chuẩn hóa cơ sở (loại trừ T21/T18/T13 để tránh nhiễu)
NORM_CHRS_BASE = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr14", "chr15", "chr16",
    "chr17", "chr19", "chr20", "chr22"
]


def compute_chr_representation(df: pd.DataFrame, target_chr: str,
                                norm_chrs: list) -> float:
    """
    EN: Compute chromosome representation with a custom normalizer list.
        The target chromosome is excluded from its own normalizer set.
    VI: Tính tỷ lệ đại diện NST với danh sách chuẩn hóa tùy chỉnh.
        NST mục tiêu được loại trừ khỏi tập chuẩn hóa của chính nó.
    """
    valid      = ~df["is_masked"] & df["norm_count"].notna()
    target_sum = df[valid & (df["chr"] == target_chr)]["norm_count"].sum()
    norm_sum   = df[valid & df["chr"].isin(norm_chrs)]["norm_count"].sum()
    if norm_sum == 0:
        return float("nan")
    return float(target_sum / norm_sum)


def robust_zscore(value: float, ref_values: np.ndarray) -> float:
    """
    EN: Robust Z-score (MAD-based). Falls back to SD when MAD=0.
    VI: Z-score bền vững (dựa trên MAD). Dùng SD khi MAD=0.
    """
    med = float(np.median(ref_values))
    mad = float(np.median(np.abs(ref_values - med)))
    if mad == 0:
        std = float(np.std(ref_values))
        return (value - med) / std if std > 0 else 0.0
    return (value - med) / (mad / 0.6745)


def main():
    # EN: Parse command-line arguments
    # VI: Phân tích đối số dòng lệnh
    parser = argparse.ArgumentParser(
        description="EN: Rare Autosomal Aneuploidy scan | VI: Quét lệch bội NST thường hiếm gặp"
    )
    parser.add_argument("--norm-counts",       required=True,
                        help="EN: Normalized bin counts Parquet | VI: File Parquet bins đã chuẩn hóa")
    parser.add_argument("--ff-metrics",        required=True,
                        help="EN: Fetal fraction JSON   | VI: File JSON tỷ lệ DNA thai nhi")
    parser.add_argument("--euploid-ref",       required=True,
                        help="EN: Euploid reference     | VI: Tham chiếu lưỡng bội")
    parser.add_argument("--sample-id",         required=True,
                        help="EN: Sample identifier     | VI: Mã mẫu")
    parser.add_argument("--chromosomes",
                        default=",".join(NORM_CHRS_BASE),
                        help="EN: Chromosomes to scan (csv) | VI: Danh sách NST cần quét (phân cách phẩy)")
    parser.add_argument("--zscore-threshold",  type=float, default=4.0,
                        help="EN: Z-score threshold (stricter than T21) | "
                             "VI: Ngưỡng Z-score (nghiêm ngặt hơn T21)")
    parser.add_argument("--min-bins",          type=int,   default=5,
                        help="EN: Minimum valid bins per chromosome | "
                             "VI: Số bins hợp lệ tối thiểu mỗi NST")
    parser.add_argument("--min-ff",            type=float, default=0.06,
                        help="EN: Min fetal fraction (stricter than tier1) | "
                             "VI: Tỷ lệ DNA thai nhi tối thiểu (nghiêm hơn tier1)")
    parser.add_argument("--out",               required=True,
                        help="EN: Output JSON file | VI: File JSON kết quả")
    args = parser.parse_args()

    # EN: Parse comma-separated chromosome list
    # VI: Tách danh sách NST từ chuỗi phân cách phẩy
    chromosomes = [c.strip() for c in args.chromosomes.split(",") if c.strip()]

    print(f"[INFO] EN: Starting RAA scan on {len(chromosomes)} chromosomes / "
          f"VI: Bắt đầu quét RAA trên {len(chromosomes)} NST "
          f"cho mẫu: {args.sample_id} ...", file=sys.stderr)

    # ── Load data / Tải dữ liệu ───────────────────────────────
    df = pq.read_table(args.norm_counts).to_pandas()
    with open(args.ff_metrics) as f:
        ff_data = json.load(f)

    ref = pq.read_table(args.euploid_ref).to_pandas() \
        if args.euploid_ref.endswith(".parquet") \
        else pd.read_csv(args.euploid_ref)

    ff_final = ff_data.get("ff_final")
    # EN: Stricter FF requirement for RAA (6% vs 4% for T21)
    # VI: Yêu cầu FF nghiêm hơn cho RAA (6% so với 4% cho T21)
    ff_ok    = ff_final is not None and ff_final >= args.min_ff

    calls = {}

    for chr_name in chromosomes:
        # EN: Exclude the target from its own normalizer set
        # VI: Loại NST mục tiêu khỏi tập chuẩn hóa của chính nó
        norm_chrs = [c for c in NORM_CHRS_BASE if c != chr_name]
        cr = compute_chr_representation(df, chr_name, norm_chrs)

        # EN: Count valid bins on this chromosome
        # VI: Đếm số bins hợp lệ trên NST này
        valid_bins = df[
            (df["chr"] == chr_name) &
            ~df["is_masked"] &
            df["norm_count"].notna()
        ]
        n_bins = len(valid_bins)

        # ── FF gate / Cổng FF ────────────────────────────────
        if not ff_ok:
            calls[chr_name] = {
                "chromosome":     chr_name,
                "cr":             round(float(cr), 6) if not np.isnan(cr) else None,
                "z":              None,
                "n_valid_bins":   n_bins,
                "classification": "NO_CALL",
                "no_call_reason": f"low_fetal_fraction:{ff_final}"
            }
            continue

        # ── Bin count gate / Cổng số bins ───────────────────
        if n_bins < args.min_bins:
            calls[chr_name] = {
                "chromosome":     chr_name,
                "cr":             round(float(cr), 6) if not np.isnan(cr) else None,
                "z":              None,
                "n_valid_bins":   n_bins,
                "classification": "NO_CALL",
                "no_call_reason": f"insufficient_bins:{n_bins}"
            }
            continue

        # ── NaN CR gate / Cổng CR NaN ────────────────────────
        if np.isnan(cr):
            calls[chr_name] = {
                "chromosome":     chr_name,
                "cr":             None,
                "z":              None,
                "n_valid_bins":   n_bins,
                "classification": "NO_CALL",
                "no_call_reason": "no_reads_in_region"
            }
            continue

        # EN: Fetch reference distribution from euploid reference matrix
        # VI: Lấy phân phối tham chiếu từ ma trận mẫu lưỡng bội
        chr_col = f"cr_{chr_name}"
        if chr_col in ref.columns:
            ref_cr = ref[chr_col].dropna().values
        elif chr_name in ref.columns:
            ref_cr = ref[chr_name].dropna().values
        else:
            ref_cr = np.array([])

        # EN: Require at least 5 reference samples for reliable Z-score
        # VI: Yêu cầu ít nhất 5 mẫu tham chiếu để Z-score đáng tin cậy
        if len(ref_cr) < 5:
            calls[chr_name] = {
                "chromosome":     chr_name,
                "cr":             round(float(cr), 6),
                "z":              None,
                "n_valid_bins":   n_bins,
                "classification": "NO_CALL",
                "no_call_reason": "insufficient_reference"
            }
            continue

        z = robust_zscore(cr, ref_cr)

        # EN: RAA flagging uses stricter |Z| ≥ 4.0 (vs 3.0 for T21)
        # VI: Đánh dấu RAA dùng |Z| ≥ 4.0 nghiêm hơn (so với 3.0 cho T21)
        classification = "SUSPECTED_RAA" if abs(z) >= args.zscore_threshold else "LOW_RISK"

        calls[chr_name] = {
            "chromosome":     chr_name,
            "cr":             round(float(cr), 6),
            "z":              round(float(z), 4),
            "n_valid_bins":   n_bins,
            "classification": classification,
            "no_call_reason": None
        }

    # ── Assemble result / Tổng hợp kết quả ───────────────────
    result = {
        "sample_id":    args.sample_id,
        "timestamp":    datetime.now(timezone.utc).isoformat(),
        "model":        "raa_robust_zscore_v1",
        "ff_final":     ff_final,
        "ff_sufficient": ff_ok,
        "calls":        calls,
        "thresholds": {
            "zscore":   args.zscore_threshold,
            "min_bins": args.min_bins,
            "min_ff":   args.min_ff
        }
    }

    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    # EN: Print summary to stderr
    # VI: In tóm tắt ra stderr
    flagged = [c for c, v in calls.items() if v["classification"] == "SUSPECTED_RAA"]
    print(f"[DONE / HOÀN THÀNH] EN: RAA scan complete / VI: Hoàn thành quét RAA "
          f"cho mẫu {args.sample_id}: "
          f"{len(flagged)} suspected / nghi ngờ RAA(s) / {len(chromosomes)} NST",
          file=sys.stderr)


if __name__ == "__main__":
    main()
