#!/usr/bin/env python3
"""
raa_scan.py
Rare Autosomal Aneuploidy (RAA) scan.

Scans specified autosomes for whole-chromosome aneuploidies
using stricter thresholds. Results are flagged with explicit
"insufficient evidence" where confidence is low.
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


NORM_CHRS_BASE = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr14", "chr15", "chr16",
    "chr17", "chr19", "chr20", "chr22"
]


def compute_chr_representation(df: pd.DataFrame, target_chr: str,
                                norm_chrs: list) -> float:
    valid      = ~df["is_masked"] & df["norm_count"].notna()
    target_sum = df[valid & (df["chr"] == target_chr)]["norm_count"].sum()
    norm_sum   = df[valid & df["chr"].isin(norm_chrs)]["norm_count"].sum()
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


def main():
    parser = argparse.ArgumentParser(description="RAA scan module")
    parser.add_argument("--norm-counts",       required=True)
    parser.add_argument("--ff-metrics",        required=True)
    parser.add_argument("--euploid-ref",       required=True)
    parser.add_argument("--sample-id",         required=True)
    parser.add_argument("--chromosomes",
                        default=",".join(NORM_CHRS_BASE))
    parser.add_argument("--zscore-threshold",  type=float, default=4.0)
    parser.add_argument("--min-bins",          type=int,   default=5)
    parser.add_argument("--min-ff",            type=float, default=0.06)
    parser.add_argument("--out",               required=True)
    args = parser.parse_args()

    chromosomes = [c.strip() for c in args.chromosomes.split(",") if c.strip()]

    df = pq.read_table(args.norm_counts).to_pandas()
    with open(args.ff_metrics) as f:
        ff_data = json.load(f)

    ref = pq.read_table(args.euploid_ref).to_pandas() \
        if args.euploid_ref.endswith(".parquet") \
        else pd.read_csv(args.euploid_ref)

    ff_final = ff_data.get("ff_final")
    ff_ok    = ff_final is not None and ff_final >= args.min_ff

    calls = {}

    for chr_name in chromosomes:
        # Use all other autosomes as normalizers (excluding target)
        norm_chrs = [c for c in NORM_CHRS_BASE if c != chr_name]
        cr = compute_chr_representation(df, chr_name, norm_chrs)

        # Check number of valid bins on this chromosome
        valid_bins = df[
            (df["chr"] == chr_name) &
            ~df["is_masked"] &
            df["norm_count"].notna()
        ]
        n_bins = len(valid_bins)

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

        # Reference distribution
        chr_col = f"cr_{chr_name}"
        if chr_col in ref.columns:
            ref_cr = ref[chr_col].dropna().values
        elif chr_name in ref.columns:
            ref_cr = ref[chr_name].dropna().values
        else:
            ref_cr = np.array([])

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

        classification = "SUSPECTED_RAA" if abs(z) >= args.zscore_threshold else "LOW_RISK"

        calls[chr_name] = {
            "chromosome":     chr_name,
            "cr":             round(float(cr), 6),
            "z":              round(float(z), 4),
            "n_valid_bins":   n_bins,
            "classification": classification,
            "no_call_reason": None
        }

    result = {
        "sample_id":   args.sample_id,
        "timestamp":   datetime.now(timezone.utc).isoformat(),
        "model":       "raa_robust_zscore_v1",
        "ff_final":    ff_final,
        "ff_sufficient": ff_ok,
        "calls":       calls,
        "thresholds": {
            "zscore":   args.zscore_threshold,
            "min_bins": args.min_bins,
            "min_ff":   args.min_ff
        }
    }

    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    flagged = [c for c, v in calls.items() if v["classification"] == "SUSPECTED_RAA"]
    print(f"RAA scan complete for {args.sample_id}: "
          f"{len(flagged)} suspected RAA(s) out of {len(chromosomes)} chromosomes",
          file=sys.stderr)


if __name__ == "__main__":
    main()
