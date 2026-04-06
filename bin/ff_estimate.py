#!/usr/bin/env python3
"""
ff_estimate.py
Fetal fraction estimation for NIPT.

Methods implemented:
1. SeqFF-inspired sex-independent method using autosomal bin profiles
2. Y-chromosome based method (for male fetuses)
3. Combined estimate with confidence interval

Output: JSON with ff_seqff, ff_y (if applicable), ff_final, ff_confidence
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
    print(f"ERROR: Missing dependency: {e}", file=sys.stderr)
    sys.exit(1)


def load_norm_counts(path: str) -> pd.DataFrame:
    """Load normalized bin counts from Parquet."""
    df = pq.read_table(path).to_pandas()
    return df


def estimate_ff_y(df: pd.DataFrame) -> Optional[float]:
    """
    Y-chromosome based fetal fraction estimate.
    FF_Y = 2 * (chrY_reads / total_reads)
    Only valid for male fetuses; returns None if chrY bins are absent/low.
    """
    y_bins  = df[df["chr"] == "chrY"]
    au_bins = df[~df["chr"].isin(["chrX", "chrY", "chrM"]) & ~df["is_masked"]]

    if y_bins.empty or au_bins.empty:
        return None

    y_sum  = y_bins["raw_count"].sum()
    au_sum = au_bins["raw_count"].sum()

    if au_sum == 0:
        return None

    # Y reads as fraction; scale by expected autosome baseline
    y_norm_bins = len(y_bins)
    au_norm_bins = len(au_bins)
    if y_norm_bins == 0:
        return None

    # Expected Y contribution for a male euploid fetus
    y_per_bin = y_sum / y_norm_bins
    au_per_bin = au_sum / au_norm_bins

    if au_per_bin == 0:
        return None

    ff_y = 2.0 * (y_per_bin / au_per_bin)
    # Clip to plausible range
    ff_y = float(np.clip(ff_y, 0.0, 1.0))
    return ff_y


def estimate_ff_seqff(df: pd.DataFrame) -> tuple[Optional[float], Optional[float]]:
    """
    SeqFF-inspired sex-independent fetal fraction estimation.

    Uses chromosome representation differences relative to expected
    diploid profiles. For clinical implementation, this should be
    trained on a validated euploid reference panel.

    Simplified implementation: uses variance of normalized bin counts
    on informative chromosomes as a proxy for FF.

    Returns: (ff_estimate, confidence_score)
    """
    # Use autosomal bins only, exclude sex chromosomes and masked bins
    auto_bins = df[
        ~df["chr"].isin(["chrX", "chrY", "chrM"]) &
        ~df["is_masked"] &
        df["norm_count"].notna()
    ].copy()

    if len(auto_bins) < 100:
        return None, None

    # Normalize counts
    median_count = auto_bins["norm_count"].median()
    if median_count <= 0:
        return None, None

    auto_bins["ratio"] = auto_bins["norm_count"] / median_count

    # Use chromosomes with relatively stable representation
    # as informative for FF estimation
    chr_stats = auto_bins.groupby("chr")["ratio"].agg(["mean", "std", "count"])
    chr_stats = chr_stats[chr_stats["count"] >= 5]

    if chr_stats.empty:
        return None, None

    # Global coefficient of variation as FF proxy
    # In real SeqFF, this uses a trained model on chr-specific features
    cv = chr_stats["std"].median() / chr_stats["mean"].median()
    if math.isnan(cv) or math.isinf(cv):
        return None, None

    # Empirical calibration (placeholder; replace with trained model coefficients)
    # CV typically correlates linearly with FF in the range 0.02-0.30
    ff_seqff = float(np.clip(cv * 2.5, 0.0, 0.60))

    # Confidence based on number of informative bins
    n_bins = len(auto_bins)
    confidence = float(np.clip((n_bins - 100) / 50000, 0.0, 1.0))

    return ff_seqff, confidence


def parse_bam_stats(bam_stats_path: str) -> dict:
    """Parse samtools flagstat output for basic metrics."""
    metrics = {
        "total_reads": 0,
        "mapped_reads": 0,
        "duplicate_reads": 0,
        "mapping_rate": 0.0,
        "dup_rate": 0.0
    }
    try:
        with open(bam_stats_path) as f:
            for line in f:
                if "in total" in line:
                    metrics["total_reads"] = int(line.split()[0])
                elif "mapped (" in line and "primary mapped" not in line:
                    metrics["mapped_reads"] = int(line.split()[0])
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
    parser = argparse.ArgumentParser(description="Fetal fraction estimation")
    parser.add_argument("--norm-counts",  required=True, help="Normalized counts Parquet")
    parser.add_argument("--bam-stats",    required=True, help="samtools flagstat output")
    parser.add_argument("--sample-id",    required=True, help="Sample ID")
    parser.add_argument("--reported-sex", default="unknown",
                        help="Reported fetal sex (male/female/unknown)")
    parser.add_argument("--out",          required=True, help="Output JSON file")
    args = parser.parse_args()

    df = load_norm_counts(args.norm_counts)
    bam_metrics = parse_bam_stats(args.bam_stats)

    # Ensure required columns exist
    for col in ["chr", "norm_count", "raw_count", "is_masked"]:
        if col not in df.columns:
            print(f"ERROR: Missing column '{col}' in norm_counts", file=sys.stderr)
            sys.exit(1)

    # Method 1: SeqFF-inspired
    ff_seqff, ff_confidence = estimate_ff_seqff(df)

    # Method 2: Y-based (applicable for male or unknown fetuses)
    ff_y = None
    if args.reported_sex in ("male", "unknown"):
        ff_y = estimate_ff_y(df)

    # Determine final FF estimate
    ff_final = None
    ff_method = "none"

    if ff_seqff is not None and ff_y is not None:
        # Both available: weighted average
        # Y-based is generally more accurate for male fetuses
        if args.reported_sex == "male":
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
        "ff_y":            round(ff_y, 4) if ff_y is not None else None,
        "ff_final":        ff_final_rounded,
        "ff_confidence":   round(ff_confidence, 4) if ff_confidence is not None else None,
        "ff_method":       ff_method,
        "reported_sex":    args.reported_sex,
        "bam_metrics":     bam_metrics
    }

    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    ff_str = f"{ff_final_rounded:.1%}" if ff_final_rounded is not None else "N/A"
    print(f"FF estimation complete for {args.sample_id}: FF={ff_str} ({ff_method})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
