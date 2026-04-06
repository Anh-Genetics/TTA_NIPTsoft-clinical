#!/usr/bin/env python3
"""
zscore_model.py
===============
EN: Aneuploidy calling for T21, T18, T13 using robust Z-score (NCV-style).
VI: Gọi lệch bội nhiễm sắc thể T21, T18, T13 bằng phương pháp Z-score bền vững (kiểu NCV).

Algorithm / Thuật toán:
1. EN: Compute chromosome representation (CR) = sum(chr_bins) / sum(reference_bins)
   VI: Tính tỷ lệ đại diện nhiễm sắc thể (CR) = tổng(bins mục tiêu) / tổng(bins chuẩn hóa)
2. EN: Compare to euploid reference distribution using robust Z-score
   VI: So sánh với phân phối tham chiếu lưỡng bội bằng Z-score bền vững
3. EN: Apply decision thresholds → HIGH_RISK / GREY_ZONE / LOW_RISK / NO_CALL
   VI: Áp dụng ngưỡng quyết định → CAO / VÙNG XÁM / THẤP / KHÔNG ĐỦ DỮ LIỆU
"""

import argparse
import json
import sys
from pathlib import Path
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


# EN: Reference chromosomes used for normalization denominator
# VI: Danh sách nhiễm sắc thể tham chiếu dùng làm mẫu số chuẩn hóa
NORM_CHRS = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr14", "chr15", "chr16",
    "chr17", "chr19", "chr20", "chr22"
]


def compute_chromosome_representation(df: pd.DataFrame, target_chr: str,
                                       norm_chrs: list) -> float:
    """
    EN: Compute chromosome representation ratio.
        CR = (reads on target chr) / (reads on normalizer chromosomes)
    VI: Tính tỷ lệ đại diện nhiễm sắc thể.
        CR = (số đọc trên NST mục tiêu) / (số đọc trên NST chuẩn hóa)
    """
    # EN: Exclude masked bins and NaN counts
    # VI: Loại bỏ các bins bị che khuất và giá trị NaN
    valid = ~df["is_masked"] & df["norm_count"].notna()
    target_reads = df[valid & (df["chr"] == target_chr)]["norm_count"].sum()
    norm_reads   = df[valid & df["chr"].isin(norm_chrs)]["norm_count"].sum()

    if norm_reads == 0:
        # EN: No usable normalizer bins → return NaN
        # VI: Không có bins chuẩn hóa hợp lệ → trả về NaN
        return float("nan")
    return float(target_reads / norm_reads)


def robust_zscore(value: float, ref_values: np.ndarray) -> float:
    """
    EN: Compute robust Z-score using median and MAD (outlier-resistant).
        Z = (x - median) / (MAD / 0.6745)
    VI: Tính Z-score bền vững dùng trung vị và MAD (kháng ngoại lệ).
        Z = (x - trung_vị) / (MAD / 0.6745)
    """
    med = float(np.median(ref_values))
    mad = float(np.median(np.abs(ref_values - med)))
    if mad == 0:
        # EN: Fall back to standard deviation when MAD is zero
        # VI: Dùng độ lệch chuẩn khi MAD bằng 0
        std = float(np.std(ref_values))
        if std == 0:
            return 0.0
        return (value - med) / std
    return (value - med) / (mad / 0.6745)


def ncv_zscore(value: float, ref_values: np.ndarray) -> float:
    """
    EN: Normalized Chromosome Value (NCV) Z-score using mean ± SD.
    VI: Z-score kiểu NCV dùng trung bình ± độ lệch chuẩn.
    """
    mu  = float(np.mean(ref_values))
    std = float(np.std(ref_values, ddof=1))
    if std == 0:
        return 0.0
    return (value - mu) / std


def classify_risk(z: float, high_thresh: float, grey_lo: float) -> str:
    """
    EN: Map Z-score to clinical classification.
    VI: Chuyển Z-score thành kết quả phân loại lâm sàng.
    """
    if z >= high_thresh:
        return "HIGH_RISK"    # EN: High risk  | VI: Nguy cơ cao
    elif z >= grey_lo:
        return "GREY_ZONE"    # EN: Borderline | VI: Vùng xám (cần theo dõi)
    else:
        return "LOW_RISK"     # EN: Low risk   | VI: Nguy cơ thấp


def load_euploid_reference(ref_path: str) -> pd.DataFrame:
    """
    EN: Load euploid reference matrix (Parquet or CSV).
    VI: Tải ma trận tham chiếu mẫu lưỡng bội (Parquet hoặc CSV).
    """
    if ref_path.endswith(".parquet"):
        return pq.read_table(ref_path).to_pandas()
    elif ref_path.endswith(".csv"):
        return pd.read_csv(ref_path)
    else:
        raise ValueError(f"Unsupported reference format / Định dạng không hỗ trợ: {ref_path}")


def main():
    # EN: Parse command-line arguments
    # VI: Phân tích đối số dòng lệnh
    parser = argparse.ArgumentParser(
        description="EN: Aneuploidy Z-score calling | VI: Gọi lệch bội bằng Z-score"
    )
    parser.add_argument("--norm-counts",        required=True,
                        help="EN: Normalized bin counts Parquet | VI: File Parquet bins đã chuẩn hóa")
    parser.add_argument("--ff-metrics",         required=True,
                        help="EN: Fetal fraction metrics JSON   | VI: File JSON tỷ lệ DNA thai nhi")
    parser.add_argument("--euploid-ref",        required=True,
                        help="EN: Euploid reference matrix      | VI: Ma trận tham chiếu lưỡng bội")
    parser.add_argument("--sample-id",          required=True,
                        help="EN: Sample identifier             | VI: Mã định danh mẫu")
    parser.add_argument("--targets",            default="chr21,chr18,chr13",
                        help="EN: Target chromosomes (csv)      | VI: NST mục tiêu (phân cách phẩy)")
    parser.add_argument("--zscore-high-risk",   type=float, default=3.0,
                        help="EN: Z-score threshold for HIGH_RISK | VI: Ngưỡng Z-score nguy cơ cao")
    parser.add_argument("--zscore-grey-zone",   type=float, default=2.5,
                        help="EN: Z-score threshold for GREY_ZONE | VI: Ngưỡng Z-score vùng xám")
    parser.add_argument("--out",                required=True,
                        help="EN: Output JSON file              | VI: File JSON kết quả đầu ra")
    args = parser.parse_args()

    # EN: Parse comma-separated target list
    # VI: Tách danh sách NST mục tiêu từ chuỗi phân cách phẩy
    targets = [t.strip() for t in args.targets.split(",")]

    # ── Load input data / Tải dữ liệu đầu vào ────────────────
    print(f"[INFO] EN: Loading normalized counts / VI: Đang tải dữ liệu bins chuẩn hóa "
          f"cho mẫu {args.sample_id} ...", file=sys.stderr)
    df  = pq.read_table(args.norm_counts).to_pandas()
    ref = load_euploid_reference(args.euploid_ref)

    with open(args.ff_metrics) as f:
        ff_data = json.load(f)

    ff_final      = ff_data.get("ff_final")
    # EN: Minimum fetal fraction hard-coded as 4%; final gate is in qc_gate
    # VI: Ngưỡng tỷ lệ DNA thai nhi tối thiểu 4%; cổng kiểm tra cuối ở qc_gate
    min_ff        = 0.04
    ff_sufficient = ff_final is not None and ff_final >= min_ff

    target_calls = {}

    for target_chr in targets:
        # ── Compute sample CR / Tính CR mẫu ──────────────────
        sample_cr = compute_chromosome_representation(df, target_chr, NORM_CHRS)

        # ── Fetch reference distribution / Lấy phân phối tham chiếu ──
        chr_col = f"cr_{target_chr}"
        if chr_col in ref.columns:
            ref_cr = ref[chr_col].dropna().values
        elif target_chr in ref.columns:
            ref_cr = ref[target_chr].dropna().values
        else:
            # EN: No matching reference column — will trigger NO_CALL
            # VI: Không tìm thấy cột tham chiếu — kết quả sẽ là NO_CALL
            ref_cr = np.array([0.01])

        # ── Guard: insufficient reference / Kiểm tra tham chiếu đủ mẫu ──
        if len(ref_cr) < 3:
            target_calls[target_chr] = {
                "chromosome":   target_chr,
                "sample_cr":    round(sample_cr, 6) if not np.isnan(sample_cr) else None,
                "z_robust":     None,
                "z_ncv":        None,
                "classification": "NO_CALL",
                "no_call_reason": "insufficient_reference_samples"
            }
            continue

        # ── Guard: no reads in target / Không có đọc trên NST mục tiêu ──
        if np.isnan(sample_cr):
            target_calls[target_chr] = {
                "chromosome":   target_chr,
                "sample_cr":    None,
                "z_robust":     None,
                "z_ncv":        None,
                "classification": "NO_CALL",
                "no_call_reason": "no_reads_in_target_region"
            }
            continue

        # ── Compute Z-scores / Tính Z-score ───────────────────
        z_robust = robust_zscore(sample_cr, ref_cr)
        z_ncv    = ncv_zscore(sample_cr, ref_cr)

        # EN: Use robust Z-score as primary; NCV as secondary confirmatory
        # VI: Dùng Z-score bền vững là chính; NCV là phụ để xác nhận
        primary_z = z_robust

        # ── Apply FF gate / Áp dụng cổng tỷ lệ DNA thai nhi ──
        if not ff_sufficient:
            classification = "NO_CALL"
            no_call_reason = f"low_fetal_fraction:{ff_final}"
        else:
            classification = classify_risk(primary_z, args.zscore_high_risk, args.zscore_grey_zone)
            no_call_reason = None

        target_calls[target_chr] = {
            "chromosome":   target_chr,
            "sample_cr":    round(float(sample_cr), 6),
            "ref_mean_cr":  round(float(np.mean(ref_cr)), 6),
            "ref_std_cr":   round(float(np.std(ref_cr, ddof=1)), 6),
            "z_robust":     round(float(z_robust), 4),
            "z_ncv":        round(float(z_ncv), 4),
            "classification": classification,
            "no_call_reason": no_call_reason
        }

    # ── Assemble result / Tổng hợp kết quả ───────────────────
    result = {
        "sample_id":   args.sample_id,
        "timestamp":   datetime.now(timezone.utc).isoformat(),
        "model":       "robust_zscore_v1",
        "ff_final":    ff_final,
        "ff_sufficient": ff_sufficient,
        "targets":     target_calls,
        "thresholds": {
            "high_risk":  args.zscore_high_risk,
            "grey_zone":  args.zscore_grey_zone
        }
    }

    # ── Write output / Ghi kết quả ────────────────────────────
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    # EN: Print per-target summary to stderr
    # VI: In tóm tắt kết quả từng NST ra stderr
    for chr_name, call in target_calls.items():
        z_str = f"{call.get('z_robust', 'N/A'):.2f}" if call.get("z_robust") is not None else "N/A"
        print(f"  [RESULT] {chr_name}: Z={z_str} => {call['classification']}", file=sys.stderr)

    print(f"[DONE / HOÀN THÀNH] Aneuploidy calling / Gọi lệch bội NST "
          f"cho mẫu: {args.sample_id}", file=sys.stderr)


if __name__ == "__main__":
    main()
