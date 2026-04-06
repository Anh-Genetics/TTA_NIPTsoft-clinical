#!/usr/bin/env python3
"""
microdeletion_call.py
=====================
EN: Microdeletion detection constrained to curated panel regions.
VI: Phát hiện vi mất đoạn NST giới hạn trong các vùng panel đã được tuyển chọn.

EN: Uses Z-score model on normalized bin counts within each panel region.
    Strictest QC requirements and conservative no-call behavior.
    A HIGH_RISK result here requires confirmation by diagnostic testing.
VI: Dùng mô hình Z-score trên số đếm bins đã chuẩn hóa trong mỗi vùng panel.
    Yêu cầu QC nghiêm ngặt nhất và hành vi thận trọng khi không gọi kết quả.
    Kết quả CAO NGUY CƠ cần được xác nhận bằng xét nghiệm chẩn đoán.

EN: Panel BED format: chr  start  end  syndrome_name  min_bins
VI: Định dạng file BED panel: chr  start  end  tên_hội_chứng  số_bins_tối_thiểu
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


def robust_zscore(value: float, ref_values: np.ndarray) -> float:
    """
    EN: Robust Z-score using MAD. Falls back to SD when MAD=0.
    VI: Z-score bền vững dùng MAD. Dùng SD khi MAD=0.
    """
    med = float(np.median(ref_values))
    mad = float(np.median(np.abs(ref_values - med)))
    if mad == 0:
        std = float(np.std(ref_values))
        return (value - med) / std if std > 0 else 0.0
    return (value - med) / (mad / 0.6745)


def load_panel_bed(panel_path: str) -> pd.DataFrame:
    """
    EN: Load curated microdeletion panel BED file.
        Expected columns: chr, start, end, syndrome, min_bins
    VI: Tải file BED panel vi mất đoạn đã được tuyển chọn.
        Các cột mong đợi: chr, start, end, tên_hội_chứng, số_bins_tối_thiểu
    """
    cols = ["chr", "start", "end", "syndrome", "min_bins"]
    try:
        df = pd.read_csv(panel_path, sep="\t", header=None,
                         names=cols, comment="#")
    except Exception:
        df = pd.read_csv(panel_path, sep="\t", comment="#")
        df.columns = cols[:len(df.columns)]
    return df


def get_region_representation(df_bins: pd.DataFrame,
                               chr_name: str, start: int, end: int) -> tuple:
    """
    EN: Compute mean normalized count for bins overlapping a panel region.
        Returns (mean_norm_count, n_valid_bins).
    VI: Tính trung bình số đếm chuẩn hóa cho các bins chồng lên vùng panel.
        Trả về (trung_bình_norm, số_bins_hợp_lệ).
    """
    region_bins = df_bins[
        (df_bins["chr"] == chr_name) &
        (df_bins["start"] >= start) &
        (df_bins["end"] <= end) &
        ~df_bins["is_masked"] &
        df_bins["norm_count"].notna()
    ]
    if region_bins.empty:
        return float("nan"), 0
    return float(region_bins["norm_count"].mean()), len(region_bins)


def main():
    # EN: Parse command-line arguments
    # VI: Phân tích đối số dòng lệnh
    parser = argparse.ArgumentParser(
        description="EN: Microdeletion calling module | VI: Module gọi vi mất đoạn NST"
    )
    parser.add_argument("--norm-counts",       required=True,
                        help="EN: Normalized bin counts Parquet | VI: File Parquet bins đã chuẩn hóa")
    parser.add_argument("--ff-metrics",        required=True,
                        help="EN: Fetal fraction JSON   | VI: File JSON tỷ lệ DNA thai nhi")
    parser.add_argument("--euploid-ref",       required=True,
                        help="EN: Euploid reference     | VI: Tham chiếu lưỡng bội")
    parser.add_argument("--panel-bed",         required=True,
                        help="EN: Microdeletion panel BED | VI: File BED panel vi mất đoạn")
    parser.add_argument("--sample-id",         required=True,
                        help="EN: Sample identifier     | VI: Mã mẫu")
    parser.add_argument("--zscore-threshold",  type=float, default=5.0,
                        help="EN: Z-score threshold (|Z| ≥ 5 for deletion) | "
                             "VI: Ngưỡng Z-score (|Z| ≥ 5 cho vi mất đoạn)")
    parser.add_argument("--min-ff",            type=float, default=0.08,
                        help="EN: Min fetal fraction (strictest: 8%) | "
                             "VI: Tỷ lệ DNA thai nhi tối thiểu (nghiêm nhất: 8%)")
    parser.add_argument("--min-bins",          type=int,   default=2,
                        help="EN: Min valid bins per panel region | "
                             "VI: Số bins hợp lệ tối thiểu mỗi vùng panel")
    parser.add_argument("--out",               required=True,
                        help="EN: Output JSON file | VI: File JSON kết quả")
    args = parser.parse_args()

    # ── Load data / Tải dữ liệu ───────────────────────────────
    df_bins = pq.read_table(args.norm_counts).to_pandas()
    with open(args.ff_metrics) as f:
        ff_data = json.load(f)

    panel = load_panel_bed(args.panel_bed)
    print(f"[INFO] EN: Loaded microdeletion panel with {len(panel)} regions / "
          f"VI: Đã tải panel vi mất đoạn gồm {len(panel)} vùng "
          f"cho mẫu: {args.sample_id} ...", file=sys.stderr)

    ref = pq.read_table(args.euploid_ref).to_pandas() \
        if args.euploid_ref.endswith(".parquet") \
        else pd.read_csv(args.euploid_ref)

    ff_final = ff_data.get("ff_final")
    # EN: Strictest FF requirement: 8% for microdeletion panel
    # VI: Yêu cầu FF nghiêm ngặt nhất: 8% cho panel vi mất đoạn
    ff_ok    = ff_final is not None and ff_final >= args.min_ff

    calls = {}

    for _, region in panel.iterrows():
        chr_name       = str(region["chr"])
        start          = int(region["start"])
        end            = int(region["end"])
        syndrome       = str(region.get("syndrome", f"{chr_name}:{start}-{end}"))
        panel_min_bins = int(region.get("min_bins", args.min_bins))
        region_key     = syndrome

        # ── FF gate / Cổng FF ────────────────────────────────
        if not ff_ok:
            calls[region_key] = {
                "syndrome":       syndrome,
                "region":         f"{chr_name}:{start}-{end}",
                "mean_norm":      None,
                "z":              None,
                "n_valid_bins":   0,
                "classification": "NO_CALL",
                "no_call_reason": f"low_fetal_fraction:{ff_final}"
            }
            continue

        # EN: Count valid bins overlapping this panel region
        # VI: Đếm bins hợp lệ nằm trong vùng panel này
        mean_norm, n_bins = get_region_representation(df_bins, chr_name, start, end)

        # ── Bins gate / Cổng bins ─────────────────────────────
        if n_bins < max(panel_min_bins, args.min_bins):
            calls[region_key] = {
                "syndrome":       syndrome,
                "region":         f"{chr_name}:{start}-{end}",
                "mean_norm":      round(mean_norm, 6) if not np.isnan(mean_norm) else None,
                "z":              None,
                "n_valid_bins":   n_bins,
                "classification": "NO_CALL",
                "no_call_reason": f"insufficient_bins:{n_bins}"
            }
            continue

        # ── Coverage gate / Cổng phủ sóng ────────────────────
        if np.isnan(mean_norm):
            calls[region_key] = {
                "syndrome":       syndrome,
                "region":         f"{chr_name}:{start}-{end}",
                "mean_norm":      None,
                "z":              None,
                "n_valid_bins":   n_bins,
                "classification": "NO_CALL",
                "no_call_reason": "no_coverage_in_region"
            }
            continue

        # EN: Look up reference distribution for this specific region
        # VI: Tra cứu phân phối tham chiếu cho vùng cụ thể này
        region_col = syndrome.replace(" ", "_").replace("/", "_")
        if region_col in ref.columns:
            ref_dist = ref[region_col].dropna().values
        else:
            # EN: No matching reference column → conservative no-call
            # VI: Không tìm thấy cột tham chiếu → thận trọng, không gọi
            ref_dist = np.array([])

        # EN: Require at least 5 reference samples for reliable Z-score
        # VI: Yêu cầu ít nhất 5 mẫu tham chiếu để Z-score đáng tin cậy
        if len(ref_dist) < 5:
            calls[region_key] = {
                "syndrome":       syndrome,
                "region":         f"{chr_name}:{start}-{end}",
                "mean_norm":      round(mean_norm, 6),
                "z":              None,
                "n_valid_bins":   n_bins,
                "classification": "NO_CALL",
                "no_call_reason": "insufficient_reference"
            }
            continue

        z = robust_zscore(mean_norm, ref_dist)

        # EN: Deletion = lower representation → negative Z; threshold |Z| ≥ 5
        # VI: Vi mất đoạn = biểu diễn thấp hơn → Z âm; ngưỡng |Z| ≥ 5
        classification = "HIGH_RISK" if z <= -args.zscore_threshold else "LOW_RISK"

        calls[region_key] = {
            "syndrome":       syndrome,
            "region":         f"{chr_name}:{start}-{end}",
            "mean_norm":      round(mean_norm, 6),
            "z":              round(float(z), 4),
            "n_valid_bins":   n_bins,
            "classification": classification,
            "no_call_reason": None
        }

    # ── Assemble result / Tổng hợp kết quả ───────────────────
    result = {
        "sample_id":    args.sample_id,
        "timestamp":    datetime.now(timezone.utc).isoformat(),
        "model":        "microdeletion_robust_zscore_v1",
        "ff_final":     ff_final,
        "ff_sufficient": ff_ok,
        "calls":        calls,
        "thresholds": {
            "zscore":   args.zscore_threshold,
            "min_ff":   args.min_ff,
            "min_bins": args.min_bins
        }
    }

    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    # EN: Print summary to stderr
    # VI: In tóm tắt ra stderr
    high_risk = [k for k, v in calls.items() if v["classification"] == "HIGH_RISK"]
    print(f"[DONE / HOÀN THÀNH] EN: Microdeletion calling complete / "
          f"VI: Hoàn thành gọi vi mất đoạn cho mẫu {args.sample_id}: "
          f"{len(high_risk)} HIGH RISK / CAO NGUY CƠ / {len(calls)} vùng panel",
          file=sys.stderr)


if __name__ == "__main__":
    main()
