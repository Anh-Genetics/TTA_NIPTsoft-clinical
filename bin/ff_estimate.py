#!/usr/bin/env python3
"""
ff_estimate.py
==============
EN: Fetal fraction estimation for NIPT.
VI: Ước tính tỷ lệ DNA thai nhi (fetal fraction) cho xét nghiệm NIPT.

EN: Methods implemented:
    1. SeqFF-inspired sex-independent method using autosomal bin profiles
    2. Y-chromosome based method (for male fetuses)
    3. Combined weighted estimate
VI: Các phương pháp đã triển khai:
    1. Phương pháp độc lập giới tính kiểu SeqFF dựa trên hồ sơ bins NST thường
    2. Phương pháp dựa trên NST Y (chỉ áp dụng cho thai nam)
    3. Ước tính kết hợp có trọng số

EN: Output: JSON with ff_seqff, ff_y (if applicable), ff_final, ff_confidence
VI: Đầu ra: JSON gồm ff_seqff, ff_y (nếu có), ff_final, ff_confidence
"""

import argparse
import json
import sys
import math
from pathlib import Path
from datetime import datetime, timezone
from typing import Optional

try:
    import pandas as pd
    import numpy as np
    import pyarrow.parquet as pq
except ImportError as e:
    # EN: Fatal — required libraries missing
    # VI: Lỗi nghiêm trọng — thiếu thư viện bắt buộc
    print(f"ERROR / LỖI: Missing dependency / Thiếu thư viện: {e}", file=sys.stderr)
    sys.exit(1)


def load_norm_counts(path: str) -> pd.DataFrame:
    """
    EN: Load normalized bin counts from Parquet.
    VI: Tải số đếm bins đã chuẩn hóa từ file Parquet.
    """
    return pq.read_table(path).to_pandas()


def estimate_ff_y(df: pd.DataFrame) -> Optional[float]:
    """
    EN: Y-chromosome based fetal fraction estimate.
        FF_Y ≈ 2 × (chrY_reads_per_bin / autosome_reads_per_bin)
        Only valid for male fetuses; returns None when chrY signal is absent.
    VI: Ước tính tỷ lệ DNA thai nhi dựa trên NST Y.
        FF_Y ≈ 2 × (số đọc chrY mỗi bin / số đọc NST thường mỗi bin)
        Chỉ hợp lệ cho thai nam; trả về None khi không có tín hiệu chrY.
    """
    y_bins  = df[df["chr"] == "chrY"]
    # EN: Autosome bins only, excluding sex chromosomes and masked bins
    # VI: Chỉ dùng bins NST thường, loại trừ NST giới tính và bins bị che khuất
    au_bins = df[~df["chr"].isin(["chrX", "chrY", "chrM"]) & ~df["is_masked"]]

    if y_bins.empty or au_bins.empty:
        return None

    y_sum  = y_bins["raw_count"].sum()
    au_sum = au_bins["raw_count"].sum()

    if au_sum == 0:
        return None

    y_norm_bins  = len(y_bins)
    au_norm_bins = len(au_bins)
    if y_norm_bins == 0:
        return None

    # EN: Normalize by bin count to get per-bin density
    # VI: Chuẩn hóa theo số lượng bin để lấy mật độ mỗi bin
    y_per_bin  = y_sum  / y_norm_bins
    au_per_bin = au_sum / au_norm_bins

    if au_per_bin == 0:
        return None

    ff_y = 2.0 * (y_per_bin / au_per_bin)
    # EN: Clip to biologically plausible range [0, 1]
    # VI: Giới hạn trong khoảng sinh học hợp lý [0, 1]
    return float(np.clip(ff_y, 0.0, 1.0))


def estimate_ff_seqff(df: pd.DataFrame) -> tuple[Optional[float], Optional[float]]:
    """
    EN: SeqFF-inspired sex-independent fetal fraction estimation.
        Uses coefficient of variation of normalized autosomal bin ratios
        as a proxy for FF. Real SeqFF requires a trained linear model;
        the calibration coefficient here is a placeholder and MUST be
        replaced with lab-validated coefficients before clinical use.
    VI: Ước tính tỷ lệ DNA thai nhi kiểu SeqFF, không phụ thuộc giới tính.
        Dùng hệ số biến thiên của tỷ lệ bins NST thường đã chuẩn hóa làm
        proxy cho FF. SeqFF thực tế cần mô hình tuyến tính đã huấn luyện;
        hệ số hiệu chuẩn dưới đây là giả định và PHẢI được thay thế bằng
        hệ số đã xác nhận bởi phòng lab trước khi sử dụng lâm sàng.

    Returns / Trả về: (ff_estimate, confidence_score)
    """
    # EN: Use autosomal bins only; exclude sex chromosomes and masked/NaN entries
    # VI: Chỉ dùng bins NST thường; loại trừ NST giới tính, bins bị che khuất và NaN
    auto_bins = df[
        ~df["chr"].isin(["chrX", "chrY", "chrM"]) &
        ~df["is_masked"] &
        df["norm_count"].notna()
    ].copy()

    if len(auto_bins) < 100:
        # EN: Insufficient bins for reliable estimation
        # VI: Không đủ bins để ước tính đáng tin cậy
        return None, None

    median_count = auto_bins["norm_count"].median()
    if median_count <= 0:
        return None, None

    # EN: Compute per-bin ratio relative to genome median
    # VI: Tính tỷ lệ mỗi bin so với trung vị toàn bộ hệ gen
    auto_bins["ratio"] = auto_bins["norm_count"] / median_count

    # EN: Compute per-chromosome statistics; require at least 5 bins per chr
    # VI: Tính thống kê cho từng NST; yêu cầu ít nhất 5 bins mỗi NST
    chr_stats = auto_bins.groupby("chr")["ratio"].agg(["mean", "std", "count"])
    chr_stats = chr_stats[chr_stats["count"] >= 5]

    if chr_stats.empty:
        return None, None

    # EN: Global CV as FF proxy; calibration coefficient is lab-specific
    # VI: CV toàn cục làm proxy FF; hệ số hiệu chuẩn phụ thuộc vào phòng lab
    cv = chr_stats["std"].median() / chr_stats["mean"].median()
    if math.isnan(cv) or math.isinf(cv):
        return None, None

    # EN: Empirical calibration — PLACEHOLDER; replace with trained model
    # VI: Hiệu chuẩn thực nghiệm — GIÁ TRỊ GIỮ CHỖ; cần thay bằng mô hình đã huấn luyện
    ff_seqff = float(np.clip(cv * 2.5, 0.0, 0.60))

    # EN: Confidence scales with number of usable autosomal bins
    # VI: Độ tin cậy tăng theo số lượng bins NST thường có thể dùng được
    n_bins = len(auto_bins)
    confidence = float(np.clip((n_bins - 100) / 50000, 0.0, 1.0))

    return ff_seqff, confidence


def parse_bam_stats(bam_stats_path: str) -> dict:
    """
    EN: Parse samtools flagstat output for basic BAM metrics.
    VI: Phân tích file đầu ra của samtools flagstat để lấy các chỉ số BAM cơ bản.
    """
    metrics = {
        "total_reads":     0,
        "mapped_reads":    0,
        "duplicate_reads": 0,
        "mapping_rate":    0.0,
        "dup_rate":        0.0
    }
    try:
        with open(bam_stats_path) as f:
            for line in f:
                if "in total" in line:
                    metrics["total_reads"]     = int(line.split()[0])
                elif "mapped (" in line and "primary mapped" not in line:
                    metrics["mapped_reads"]    = int(line.split()[0])
                elif "duplicates" in line:
                    metrics["duplicate_reads"] = int(line.split()[0])
        if metrics["total_reads"] > 0:
            metrics["mapping_rate"] = metrics["mapped_reads"] / metrics["total_reads"]
        if metrics["mapped_reads"] > 0:
            metrics["dup_rate"] = metrics["duplicate_reads"] / metrics["mapped_reads"]
    except (FileNotFoundError, ValueError, IndexError):
        pass
    return metrics


def main():
    # EN: Parse command-line arguments
    # VI: Phân tích đối số dòng lệnh
    parser = argparse.ArgumentParser(
        description="EN: Fetal fraction estimation | VI: Ước tính tỷ lệ DNA thai nhi"
    )
    parser.add_argument("--norm-counts",  required=True,
                        help="EN: Normalized counts Parquet | VI: File Parquet bins đã chuẩn hóa")
    parser.add_argument("--bam-stats",    required=True,
                        help="EN: samtools flagstat output  | VI: File đầu ra flagstat")
    parser.add_argument("--sample-id",    required=True,
                        help="EN: Sample identifier         | VI: Mã định danh mẫu")
    parser.add_argument("--reported-sex", default="unknown",
                        help="EN: Reported fetal sex (male/female/unknown) | "
                             "VI: Giới tính thai nhi được khai báo")
    parser.add_argument("--out",          required=True,
                        help="EN: Output JSON file          | VI: File JSON kết quả đầu ra")
    args = parser.parse_args()

    print(f"[INFO] EN: Starting fetal fraction estimation / "
          f"VI: Bắt đầu ước tính tỷ lệ DNA thai nhi cho mẫu: {args.sample_id} ...",
          file=sys.stderr)

    df = load_norm_counts(args.norm_counts)
    bam_metrics = parse_bam_stats(args.bam_stats)

    # EN: Verify required columns exist
    # VI: Kiểm tra các cột bắt buộc tồn tại trong dữ liệu
    for col in ["chr", "norm_count", "raw_count", "is_masked"]:
        if col not in df.columns:
            print(f"ERROR / LỖI: Missing column / Thiếu cột '{col}' in norm_counts",
                  file=sys.stderr)
            sys.exit(1)

    # ── Method 1: SeqFF-inspired / Phương pháp 1: Kiểu SeqFF ─
    print(f"[INFO] EN: Running SeqFF-style estimation / "
          f"VI: Chạy ước tính kiểu SeqFF ...", file=sys.stderr)
    ff_seqff, ff_confidence = estimate_ff_seqff(df)

    # ── Method 2: Y-based (male/unknown) / Phương pháp 2: Dựa trên chrY ─
    ff_y = None
    if args.reported_sex in ("male", "unknown"):
        print(f"[INFO] EN: Running Y-chromosome FF estimation / "
              f"VI: Chạy ước tính FF dựa trên NST Y ...", file=sys.stderr)
        ff_y = estimate_ff_y(df)

    # ── Combine estimates / Kết hợp các ước tính ─────────────
    ff_final = None
    ff_method = "none"

    if ff_seqff is not None and ff_y is not None:
        # EN: Both available: weighted by reported sex
        # VI: Cả hai phương pháp có kết quả: kết hợp có trọng số theo giới tính
        if args.reported_sex == "male":
            # EN: Y-based is generally more accurate for confirmed male fetuses
            # VI: Phương pháp Y chính xác hơn cho thai nam đã xác nhận
            ff_final = 0.4 * ff_seqff + 0.6 * ff_y
        else:
            ff_final = 0.7 * ff_seqff + 0.3 * ff_y
        ff_method = "combined"
    elif ff_seqff is not None:
        ff_final = ff_seqff
        ff_method = "seqff"
    elif ff_y is not None:
        ff_final = ff_y
        ff_method = "y_based"

    ff_final_rounded = round(ff_final, 4) if ff_final is not None else None

    result = {
        "sample_id":       args.sample_id,
        "timestamp":       datetime.now(timezone.utc).isoformat(),
        "ff_seqff":        round(ff_seqff, 4) if ff_seqff is not None else None,
        "ff_y":            round(ff_y, 4)     if ff_y     is not None else None,
        "ff_final":        ff_final_rounded,
        "ff_confidence":   round(ff_confidence, 4) if ff_confidence is not None else None,
        "ff_method":       ff_method,
        "reported_sex":    args.reported_sex,
        "bam_metrics":     bam_metrics
    }

    # ── Write output / Ghi kết quả ────────────────────────────
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    ff_str = f"{ff_final_rounded:.1%}" if ff_final_rounded is not None else "N/A"
    print(f"[DONE / HOÀN THÀNH] EN: Fetal fraction estimation complete / "
          f"VI: Hoàn thành ước tính tỷ lệ DNA thai nhi "
          f"cho mẫu {args.sample_id}: FF={ff_str} ({ff_method})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
