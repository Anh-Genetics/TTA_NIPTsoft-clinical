#!/usr/bin/env python3
"""
microdeletion_call.py
Microdeletion detection constrained to curated panel regions.

Uses Z-score model on normalized bin counts within each panel region.
Strictest QC requirements and conservative no-call behavior.

Panel BED format: chr  start  end  syndrome_name  min_bins
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


def robust_zscore(value: float, ref_values: np.ndarray) -> float:
    med = float(np.median(ref_values))
    mad = float(np.median(np.abs(ref_values - med)))
    if mad == 0:
        std = float(np.std(ref_values))
        return (value - med) / std if std > 0 else 0.0
    return (value - med) / (mad / 0.6745)


def load_panel_bed(panel_path: str) -> pd.DataFrame:
    """Load curated microdeletion panel BED."""
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
    Compute mean normalized count for bins overlapping a region.
    Returns (mean_norm, n_valid_bins).
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
    parser = argparse.ArgumentParser(description="Microdeletion calling module")
    parser.add_argument("--norm-counts",       required=True)
    parser.add_argument("--ff-metrics",        required=True)
    parser.add_argument("--euploid-ref",       required=True)
    parser.add_argument("--panel-bed",         required=True)
    parser.add_argument("--sample-id",         required=True)
    parser.add_argument("--zscore-threshold",  type=float, default=5.0)
    parser.add_argument("--min-ff",            type=float, default=0.08)
    parser.add_argument("--min-bins",          type=int,   default=2)
    parser.add_argument("--out",               required=True)
    args = parser.parse_args()

    df_bins = pq.read_table(args.norm_counts).to_pandas()
    with open(args.ff_metrics) as f:
        ff_data = json.load(f)

    panel = load_panel_bed(args.panel_bed)

    ref = pq.read_table(args.euploid_ref).to_pandas() \
        if args.euploid_ref.endswith(".parquet") \
        else pd.read_csv(args.euploid_ref)

    ff_final = ff_data.get("ff_final")
    ff_ok    = ff_final is not None and ff_final >= args.min_ff

    calls = {}

    for _, region in panel.iterrows():
        chr_name     = str(region["chr"])
        start        = int(region["start"])
        end          = int(region["end"])
        syndrome     = str(region.get("syndrome", f"{chr_name}:{start}-{end}"))
        panel_min_bins = int(region.get("min_bins", args.min_bins))

        region_key = syndrome

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

        mean_norm, n_bins = get_region_representation(df_bins, chr_name, start, end)

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

        # Reference distribution for this region
        region_col = syndrome.replace(" ", "_").replace("/", "_")
        if region_col in ref.columns:
            ref_dist = ref[region_col].dropna().values
        else:
            # Fallback: use genome-wide mean
            ref_dist = np.array([])

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

        # Deletion = lower representation => negative Z
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

    high_risk = [k for k, v in calls.items() if v["classification"] == "HIGH_RISK"]
    print(f"Microdeletion calling complete for {args.sample_id}: "
          f"{len(high_risk)} HIGH RISK region(s) out of {len(calls)} panel regions",
          file=sys.stderr)


if __name__ == "__main__":
    main()
