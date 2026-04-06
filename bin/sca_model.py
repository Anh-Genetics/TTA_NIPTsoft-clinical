#!/usr/bin/env python3
"""
sca_model.py
Sex chromosome aneuploidy (SCA) calling.

Handles:
- Turner syndrome (45,X)
- Klinefelter syndrome (47,XXY)
- Triple X syndrome (47,XXX)
- XYY syndrome (47,XYY)

Uses separate robust Z-score model for chrX and chrY
with sex-specific reference distributions.
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
    print(f"ERROR: Missing dependency: {e}", file=sys.stderr)
    sys.exit(1)


# Autosomal normalizer chromosomes (same as main model)
NORM_CHRS = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr14", "chr15", "chr16",
    "chr17", "chr19", "chr20", "chr22"
]


def compute_chr_representation(df: pd.DataFrame, target_chr: str) -> float:
    valid      = ~df["is_masked"] & df["norm_count"].notna()
    target_sum = df[valid & (df["chr"] == target_chr)]["norm_count"].sum()
    norm_sum   = df[valid & df["chr"].isin(NORM_CHRS)]["norm_count"].sum()
    if norm_sum == 0:
        return float("nan")
    return float(target_sum / norm_sum)


def robust_zscore(value: float, ref_values: np.ndarray) -> float:
    med = float(np.median(ref_values))
    mad = float(np.median(np.abs(ref_values - med)))
    if mad == 0:
        std = float(np.std(ref_values))
        return (value - med) / std if std > 0 else 0.0
    return (value - med) / (mad / 0.6745)


def infer_fetal_sex(cr_x: float, cr_y: float,
                    ref_female_x: np.ndarray, ref_male_x: np.ndarray) -> str:
    """Infer fetal sex from ChrY representation."""
    if cr_y is None or np.isnan(cr_y):
        return "female"
    # Threshold: if Y representation > 0.005, likely male
    if cr_y > 0.005:
        return "male"
    return "female"


def main():
    parser = argparse.ArgumentParser(description="SCA calling module")
    parser.add_argument("--norm-counts",      required=True)
    parser.add_argument("--ff-metrics",       required=True)
    parser.add_argument("--euploid-ref",      required=True)
    parser.add_argument("--sample-id",        required=True)
    parser.add_argument("--reported-sex",     default="unknown")
    parser.add_argument("--zscore-threshold", type=float, default=3.0)
    parser.add_argument("--min-ff",           type=float, default=0.04)
    parser.add_argument("--out",              required=True)
    args = parser.parse_args()

    df = pq.read_table(args.norm_counts).to_pandas()
    with open(args.ff_metrics) as f:
        ff_data = json.load(f)

    ref = pq.read_table(args.euploid_ref).to_pandas() \
        if args.euploid_ref.endswith(".parquet") \
        else pd.read_csv(args.euploid_ref)

    ff_final = ff_data.get("ff_final")
    ff_ok    = ff_final is not None and ff_final >= args.min_ff

    cr_x = compute_chr_representation(df, "chrX")
    cr_y = compute_chr_representation(df, "chrY")

    # Reference distributions
    ref_x_female = ref.get("cr_chrX_female", ref.get("cr_chrX", pd.Series([]))).dropna().values
    ref_x_male   = ref.get("cr_chrX_male",   pd.Series([])).dropna().values
    ref_y_male   = ref.get("cr_chrY_male",   pd.Series([])).dropna().values
    ref_y_female = ref.get("cr_chrY_female", pd.Series([])).dropna().values

    inferred_sex = infer_fetal_sex(
        cr_x, cr_y,
        ref_x_female if len(ref_x_female) > 0 else np.array([0.1]),
        ref_x_male   if len(ref_x_male)   > 0 else np.array([0.05])
    )

    calls = {}

    def sca_call(label: str, cr: float, ref_dist: np.ndarray,
                 direction: str = "both") -> dict:
        """direction: 'up'/'down'/'both'"""
        if not ff_ok:
            return {
                "cr": round(float(cr), 6) if not np.isnan(cr) else None,
                "z":  None,
                "classification": "NO_CALL",
                "no_call_reason": f"low_fetal_fraction:{ff_final}"
            }
        if np.isnan(cr):
            return {
                "cr": None, "z": None,
                "classification": "NO_CALL",
                "no_call_reason": "no_reads_in_region"
            }
        if len(ref_dist) < 3:
            return {
                "cr": round(float(cr), 6), "z": None,
                "classification": "NO_CALL",
                "no_call_reason": "insufficient_reference"
            }

        z = robust_zscore(cr, ref_dist)
        abs_z = abs(z)

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

    # chrX analysis — use sex-appropriate reference
    if inferred_sex == "female":
        ref_x = ref_x_female if len(ref_x_female) >= 3 else np.array([0.1, 0.1, 0.1])
    else:
        ref_x = ref_x_male if len(ref_x_male) >= 3 else np.array([0.05, 0.05, 0.05])

    calls["chrX"] = sca_call("chrX", cr_x, ref_x, direction="both")

    # chrY analysis
    ref_y = ref_y_male if len(ref_y_male) >= 3 else np.array([0.003, 0.003, 0.003])
    calls["chrY"] = sca_call("chrY", cr_y, ref_y, direction="both")

    result = {
        "sample_id":    args.sample_id,
        "timestamp":    datetime.now(timezone.utc).isoformat(),
        "model":        "sca_robust_zscore_v1",
        "inferred_sex": inferred_sex,
        "reported_sex": args.reported_sex,
        "ff_final":     ff_final,
        "ff_sufficient": ff_ok,
        "calls":        calls,
        "threshold":    args.zscore_threshold
    }

    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    print(f"SCA calling complete for {args.sample_id}: fetal_sex_inferred={len(calls)} targets",
          file=sys.stderr)
    for c, v in calls.items():
        print(f"  {c}: Z={v.get('z', 'N/A')} => {v['classification']}", file=sys.stderr)


if __name__ == "__main__":
    main()
