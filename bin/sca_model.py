#!/usr/bin/env python3
"""
sca_model.py
============
EN: Sex chromosome aneuploidy (SCA) calling module.
VI: Module gọi lệch bội nhiễm sắc thể giới tính (SCA).

EN: Handles:
    - Turner syndrome  (45,X)  — low chrX, absent chrY
    - Klinefelter syndrome (47,XXY) — normal/high chrX, chrY present
    - Triple X syndrome (47,XXX) — elevated chrX
    - XYY syndrome (47,XYY)    — elevated chrY
VI: Xử lý:
    - Hội chứng Turner  (45,X)  — chrX thấp, không có chrY
    - Hội chứng Klinefelter (47,XXY) — chrX bình thường/cao, có chrY
    - Hội chứng Triple X (47,XXX) — chrX tăng cao
    - Hội chứng XYY (47,XYY)    — chrY tăng cao

EN: Uses separate robust Z-score model for chrX and chrY
    with sex-specific reference distributions.
VI: Dùng mô hình Z-score bền vững riêng cho chrX và chrY
    với phân phối tham chiếu riêng theo giới tính.
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


# EN: Autosomal normalizer chromosomes (same list as main aneuploidy model)
# VI: Danh sách NST thường dùng để chuẩn hóa (giống mô hình lệch bội chính)
NORM_CHRS = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr14", "chr15", "chr16",
    "chr17", "chr19", "chr20", "chr22"
]


def compute_chr_representation(df: pd.DataFrame, target_chr: str) -> float:
    """
    EN: Compute representation of a sex chromosome relative to autosomes.
    VI: Tính tỷ lệ đại diện của NST giới tính so với NST thường.
    """
    valid      = ~df["is_masked"] & df["norm_count"].notna()
    target_sum = df[valid & (df["chr"] == target_chr)]["norm_count"].sum()
    norm_sum   = df[valid & df["chr"].isin(NORM_CHRS)]["norm_count"].sum()
    if norm_sum == 0:
        return float("nan")
    return float(target_sum / norm_sum)


def robust_zscore(value: float, ref_values: np.ndarray) -> float:
    """
    EN: Robust Z-score using MAD (median absolute deviation).
    VI: Z-score bền vững dùng MAD (độ lệch tuyệt đối trung vị).
    """
    med = float(np.median(ref_values))
    mad = float(np.median(np.abs(ref_values - med)))
    if mad == 0:
        std = float(np.std(ref_values))
        return (value - med) / std if std > 0 else 0.0
    return (value - med) / (mad / 0.6745)


def infer_fetal_sex(cr_x: float, cr_y: float,
                    ref_female_x: np.ndarray, ref_male_x: np.ndarray) -> str:
    """
    EN: Infer fetal sex from chrY representation.
        Threshold: cr_y > 0.005 → likely male fetus.
    VI: Suy luận giới tính thai nhi từ tỷ lệ đại diện chrY.
        Ngưỡng: cr_y > 0.005 → khả năng thai nam.
    """
    if cr_y is None or np.isnan(cr_y):
        # EN: No chrY signal → default to female
        # VI: Không có tín hiệu chrY → mặc định là thai nữ
        return "female"
    if cr_y > 0.005:
        return "male"
    return "female"


def main():
    # EN: Parse command-line arguments
    # VI: Phân tích đối số dòng lệnh
    parser = argparse.ArgumentParser(
        description="EN: SCA calling module | VI: Module gọi lệch bội NST giới tính"
    )
    parser.add_argument("--norm-counts",      required=True,
                        help="EN: Normalized bin counts Parquet | VI: File Parquet bins đã chuẩn hóa")
    parser.add_argument("--ff-metrics",       required=True,
                        help="EN: Fetal fraction JSON   | VI: File JSON tỷ lệ DNA thai nhi")
    parser.add_argument("--euploid-ref",      required=True,
                        help="EN: Euploid reference     | VI: Tham chiếu lưỡng bội")
    parser.add_argument("--sample-id",        required=True,
                        help="EN: Sample identifier     | VI: Mã mẫu")
    parser.add_argument("--reported-sex",     default="unknown",
                        help="EN: Reported fetal sex    | VI: Giới tính thai nhi được khai báo")
    parser.add_argument("--zscore-threshold", type=float, default=3.0,
                        help="EN: Z-score threshold     | VI: Ngưỡng Z-score")
    parser.add_argument("--min-ff",           type=float, default=0.04,
                        help="EN: Min fetal fraction    | VI: Tỷ lệ DNA thai nhi tối thiểu")
    parser.add_argument("--out",              required=True,
                        help="EN: Output JSON file      | VI: File JSON kết quả")
    args = parser.parse_args()

    print(f"[INFO] EN: Starting SCA calling / VI: Bắt đầu gọi lệch bội NST giới tính "
          f"cho mẫu: {args.sample_id} ...", file=sys.stderr)

    # ── Load data / Tải dữ liệu ───────────────────────────────
    df = pq.read_table(args.norm_counts).to_pandas()
    with open(args.ff_metrics) as f:
        ff_data = json.load(f)

    ref = pq.read_table(args.euploid_ref).to_pandas() \
        if args.euploid_ref.endswith(".parquet") \
        else pd.read_csv(args.euploid_ref)

    ff_final = ff_data.get("ff_final")
    ff_ok    = ff_final is not None and ff_final >= args.min_ff

    # ── Compute sex chromosome representations / Tính CR NST giới tính ──
    cr_x = compute_chr_representation(df, "chrX")
    cr_y = compute_chr_representation(df, "chrY")

    # EN: Retrieve sex-stratified reference distributions
    # VI: Lấy phân phối tham chiếu phân tầng theo giới tính
    ref_x_female = ref.get("cr_chrX_female", ref.get("cr_chrX", pd.Series([]))).dropna().values
    ref_x_male   = ref.get("cr_chrX_male",   pd.Series([])).dropna().values
    ref_y_male   = ref.get("cr_chrY_male",   pd.Series([])).dropna().values
    ref_y_female = ref.get("cr_chrY_female", pd.Series([])).dropna().values

    # EN: Infer fetal sex to choose the correct reference
    # VI: Suy luận giới tính thai nhi để chọn tham chiếu phù hợp
    inferred_sex = infer_fetal_sex(
        cr_x, cr_y,
        ref_x_female if len(ref_x_female) > 0 else np.array([0.1]),
        ref_x_male   if len(ref_x_male)   > 0 else np.array([0.05])
    )

    calls = {}

    def sca_call(label: str, cr: float, ref_dist: np.ndarray,
                 direction: str = "both") -> dict:
        """
        EN: Compute a single SCA call for one chromosome.
            direction: 'up' (gain), 'down' (loss), 'both'
        VI: Tính kết quả SCA cho một NST.
            direction: 'up' (tăng số lượng), 'down' (mất), 'both' (cả hai)
        """
        # EN: FF gate — no call if fetal fraction is insufficient
        # VI: Cổng FF — không gọi nếu tỷ lệ DNA thai nhi không đủ
        if not ff_ok:
            return {
                "cr": round(float(cr), 6) if not np.isnan(cr) else None,
                "z":  None,
                "classification": "NO_CALL",
                "no_call_reason": f"low_fetal_fraction:{ff_final}"
            }
        # EN: No reads in region → no call
        # VI: Không có đọc trong vùng → không gọi
        if np.isnan(cr):
            return {
                "cr": None, "z": None,
                "classification": "NO_CALL",
                "no_call_reason": "no_reads_in_region"
            }
        # EN: Reference too small → no call
        # VI: Tham chiếu quá nhỏ → không gọi
        if len(ref_dist) < 3:
            return {
                "cr": round(float(cr), 6), "z": None,
                "classification": "NO_CALL",
                "no_call_reason": "insufficient_reference"
            }

        z     = robust_zscore(cr, ref_dist)
        abs_z = abs(z)

        # EN: Apply directional threshold for sex chromosomes
        # VI: Áp dụng ngưỡng có hướng cho NST giới tính
        if direction == "down":
            significant = z <= -args.zscore_threshold
        elif direction == "up":
            significant = z >= args.zscore_threshold
        else:
            significant = abs_z >= args.zscore_threshold

        return {
            "cr":             round(float(cr), 6),
            "z":              round(float(z), 4),
            "classification": "HIGH_RISK" if significant else "LOW_RISK",
            "no_call_reason": None
        }

    # ── chrX analysis / Phân tích chrX ───────────────────────
    # EN: Use sex-appropriate reference distribution for chrX
    # VI: Dùng phân phối tham chiếu phù hợp giới tính cho chrX
    if inferred_sex == "female":
        ref_x = ref_x_female if len(ref_x_female) >= 3 else np.array([0.1, 0.1, 0.1])
    else:
        ref_x = ref_x_male   if len(ref_x_male)   >= 3 else np.array([0.05, 0.05, 0.05])

    calls["chrX"] = sca_call("chrX", cr_x, ref_x, direction="both")

    # ── chrY analysis / Phân tích chrY ───────────────────────
    ref_y = ref_y_male if len(ref_y_male) >= 3 else np.array([0.003, 0.003, 0.003])
    calls["chrY"] = sca_call("chrY", cr_y, ref_y, direction="both")

    # ── Assemble result / Tổng hợp kết quả ───────────────────
    result = {
        "sample_id":     args.sample_id,
        "timestamp":     datetime.now(timezone.utc).isoformat(),
        "model":         "sca_robust_zscore_v1",
        "inferred_sex":  inferred_sex,
        "reported_sex":  args.reported_sex,
        "ff_final":      ff_final,
        "ff_sufficient": ff_ok,
        "calls":         calls,
        "threshold":     args.zscore_threshold
    }

    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    # EN: Print per-chromosome summary to stderr
    # VI: In tóm tắt kết quả từng NST ra stderr
    print(f"[INFO] EN: SCA calling complete for {len(calls)} targets / "
          f"VI: Hoàn thành gọi SCA cho {len(calls)} NST mục tiêu", file=sys.stderr)
    for c, v in calls.items():
        print(f"  [RESULT] {c}: Z={v.get('z', 'N/A')} => {v['classification']}",
              file=sys.stderr)

    print(f"[DONE / HOÀN THÀNH] SCA calling / Gọi SCA "
          f"cho mẫu: {args.sample_id}", file=sys.stderr)


if __name__ == "__main__":
    main()
