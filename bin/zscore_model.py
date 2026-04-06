#!/usr/bin/env python3
"""
zscore_model.py
Aneuploidy calling for T21, T18, T13 using robust Z-score (NCV-style).

Algorithm:
1. Compute chromosome representation (CR) = sum(chr_bins) / sum(reference_bins)
2. Compare to euploid reference distribution using robust Z-score
3. Apply decision thresholds to classify: HIGH_RISK / GREY_ZONE / LOW_RISK / NO_CALL
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
    print(f"ERROR: Missing dependency: {e}", file=sys.stderr)
    sys.exit(1)


# Reference chromosomes used for normalization denominator
NORM_CHRS = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr14", "chr15", "chr16",
    "chr17", "chr19", "chr20", "chr22"
]


def compute_chromosome_representation(df: pd.DataFrame, target_chr: str,
                                       norm_chrs: list) -> float:
    """
    Compute chromosome representation ratio.
    CR = (reads on target chr) / (reads on normalizer chromosomes)
    """
    valid = ~df["is_masked"] & df["norm_count"].notna()
    target_reads = df[valid & (df["chr"] == target_chr)]["norm_count"].sum()
    norm_reads   = df[valid & df["chr"].isin(norm_chrs)]["norm_count"].sum()

    if norm_reads == 0:
        return float("nan")
    return float(target_reads / norm_reads)


def robust_zscore(value: float, ref_values: np.ndarray) -> float:
    """
    Compute robust Z-score using median and MAD.
    Z = (x - median) / (MAD / 0.6745)
    """
    med = float(np.median(ref_values))
    mad = float(np.median(np.abs(ref_values - med)))
    if mad == 0:
        std = float(np.std(ref_values))
        if std == 0:
            return 0.0
        return (value - med) / std
    return (value - med) / (mad / 0.6745)


def ncv_zscore(value: float, ref_values: np.ndarray) -> float:
    """
    Normalized Chromosome Value (NCV) Z-score using mean ± SD.
    """
    mu  = float(np.mean(ref_values))
    std = float(np.std(ref_values, ddof=1))
    if std == 0:
        return 0.0
    return (value - mu) / std


def classify_risk(z: float, high_thresh: float, grey_lo: float) -> str:
    """Map Z-score to clinical classification."""
    if z >= high_thresh:
        return "HIGH_RISK"
    elif z >= grey_lo:
        return "GREY_ZONE"
    else:
        return "LOW_RISK"


def load_euploid_reference(ref_path: str) -> pd.DataFrame:
    """Load euploid reference matrix."""
    if ref_path.endswith(".parquet"):
        return pq.read_table(ref_path).to_pandas()
    elif ref_path.endswith(".csv"):
        return pd.read_csv(ref_path)
    else:
        raise ValueError(f"Unsupported reference format: {ref_path}")


def main():
    parser = argparse.ArgumentParser(description="Aneuploidy Z-score calling")
    parser.add_argument("--norm-counts",        required=True)
    parser.add_argument("--ff-metrics",         required=True)
    parser.add_argument("--euploid-ref",        required=True)
    parser.add_argument("--sample-id",          required=True)
    parser.add_argument("--targets",            default="chr21,chr18,chr13")
    parser.add_argument("--zscore-high-risk",   type=float, default=3.0)
    parser.add_argument("--zscore-grey-zone",   type=float, default=2.5)
    parser.add_argument("--out",                required=True)
    args = parser.parse_args()

    targets = [t.strip() for t in args.targets.split(",")]

    # Load data
    df  = pq.read_table(args.norm_counts).to_pandas()
    ref = load_euploid_reference(args.euploid_ref)

    with open(args.ff_metrics) as f:
        ff_data = json.load(f)

    ff_final      = ff_data.get("ff_final")
    min_ff        = 0.04  # hardcoded minimum; param in qc_gate
    ff_sufficient = ff_final is not None and ff_final >= min_ff

    target_calls = {}

    for target_chr in targets:
        # Sample chromosome representation
        sample_cr = compute_chromosome_representation(df, target_chr, NORM_CHRS)

        # Reference distribution for this chromosome
        chr_col = f"cr_{target_chr}"
        if chr_col in ref.columns:
            ref_cr = ref[chr_col].dropna().values
        elif target_chr in ref.columns:
            ref_cr = ref[target_chr].dropna().values
        else:
            # Fallback: compute from reference samples if they have bin data
            ref_cr = np.array([0.01])  # stub; will trigger NO_CALL

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

        z_robust = robust_zscore(sample_cr, ref_cr)
        z_ncv    = ncv_zscore(sample_cr, ref_cr)

        # Use robust Z-score as primary; NCV as secondary
        primary_z = z_robust

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

    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    for chr_name, call in target_calls.items():
        z_str = f"{call.get('z_robust', 'N/A'):.2f}" if call.get("z_robust") is not None else "N/A"
        print(f"  {chr_name}: Z={z_str} => {call['classification']}", file=sys.stderr)

    print(f"Aneuploidy calling complete for {args.sample_id}", file=sys.stderr)


if __name__ == "__main__":
    main()
