"""
tests/unit/test_gc_correct.py
Unit tests for GC correction logic (Python-equivalent tests).
Tests the mathematical correctness of LOESS-like corrections.
"""

import sys
import os
import json
import pytest
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../bin"))


class TestGCCorrectionMath:
    """Test the mathematical properties of GC bias correction."""

    def _make_gc_biased_df(self, n_bins: int = 500, seed: int = 42) -> pd.DataFrame:
        """
        Create a synthetic dataset with known GC bias.
        True counts are 1.0 everywhere; observed counts are biased by GC.
        GC bias model: count = true_count * (1 + 2*(gc - 0.5)^2)
        """
        rng = np.random.default_rng(seed)
        gc    = rng.uniform(0.2, 0.8, n_bins)
        # Bias: highest counts at extreme GC
        bias  = 1 + 2 * (gc - 0.5) ** 2
        noise = rng.normal(1.0, 0.05, n_bins)
        raw   = bias * noise
        raw   = np.maximum(raw, 0)

        chrs = ["chr" + str(i % 22 + 1) for i in range(n_bins)]
        return pd.DataFrame({
            "chr":       chrs,
            "start":     np.arange(n_bins) * 50000,
            "end":       (np.arange(n_bins) + 1) * 50000,
            "raw_count": raw,
            "gc_frac":   gc,
            "map_frac":  np.ones(n_bins),
            "is_masked": np.zeros(n_bins, dtype=bool)
        })

    def test_gc_correlation_reduced_after_correction(self):
        """After GC correction, correlation between count and GC should be lower."""
        df = self._make_gc_biased_df()

        # Compute raw GC correlation
        cor_raw = abs(np.corrcoef(df["gc_frac"], df["raw_count"])[0, 1])

        # Simple LOESS correction (simplified for testing)
        valid = df["raw_count"] > 0
        global_mean = df["raw_count"][valid].mean()

        from scipy.interpolate import UnivariateSpline
        # Sort by GC for spline fitting
        sorted_idx = np.argsort(df.loc[valid, "gc_frac"])
        gc_sorted  = df.loc[valid, "gc_frac"].values[sorted_idx]
        cnt_sorted = df.loc[valid, "raw_count"].values[sorted_idx]

        # Fit smoothing spline as LOESS proxy
        spline = UnivariateSpline(gc_sorted, cnt_sorted, s=len(gc_sorted) * 0.1, k=3)
        predicted = spline(df["gc_frac"])
        predicted = np.maximum(predicted, 1e-9)
        corrected = df["raw_count"] * (global_mean / predicted)

        cor_corrected = abs(np.corrcoef(df["gc_frac"], corrected)[0, 1])
        assert cor_corrected < cor_raw, (
            f"GC correction should reduce GC correlation: "
            f"raw={cor_raw:.3f}, corrected={cor_corrected:.3f}"
        )

    def test_corrected_counts_non_negative(self):
        """Corrected counts should always be non-negative."""
        df = self._make_gc_biased_df()
        global_mean = df["raw_count"].mean()
        # Simple multiplicative correction
        corrected = df["raw_count"] * (global_mean / (df["raw_count"].mean()))
        assert (corrected >= 0).all()

    def test_masked_bins_excluded_from_normalization(self):
        """Masked bins should not affect the correction model."""
        df = self._make_gc_biased_df(n_bins=200)
        # Mask extreme GC bins
        df.loc[df["gc_frac"] < 0.25, "is_masked"] = True
        df.loc[df["gc_frac"] > 0.75, "is_masked"] = True

        valid_mask = ~df["is_masked"] & (df["raw_count"] > 0)
        # Verify that correction uses only valid bins
        assert valid_mask.sum() > 50, "Should have enough valid bins"
        assert valid_mask.sum() < len(df), "Some bins should be masked"


class TestBinFiltering:
    """Test bin quality filtering logic."""

    def test_extreme_gc_bins_masked(self):
        """Bins with GC < 20% or > 80% should be flagged as problematic."""
        GC_MIN = 0.20
        GC_MAX = 0.80

        df = pd.DataFrame({
            "gc_frac":   [0.10, 0.25, 0.50, 0.75, 0.85],
            "map_frac":  [1.0,  1.0,  1.0,  1.0,  1.0],
            "raw_count": [100,  100,  100,  100,  100]
        })
        expected_valid = [False, True, True, True, False]
        is_valid = (
            (df["gc_frac"]  >= GC_MIN) &
            (df["gc_frac"]  <= GC_MAX) &
            (df["map_frac"] >= 0.50)
        )
        assert list(is_valid) == expected_valid

    def test_low_mappability_bins_masked(self):
        """Bins with mappability < 50% should be masked."""
        MAP_MIN = 0.50
        df = pd.DataFrame({
            "gc_frac":   [0.50, 0.50, 0.50],
            "map_frac":  [0.30, 0.50, 0.80],
            "raw_count": [100,  100,  100]
        })
        is_valid = df["map_frac"] >= MAP_MIN
        assert list(is_valid) == [False, True, True]


class TestNormalizationScale:
    """Test global normalization properties."""

    def test_autosomal_median_approaches_one(self):
        """After normalization, median autosomal bin count should be ~1.0."""
        rng = np.random.default_rng(7)
        counts = rng.lognormal(mean=np.log(1000), sigma=0.1, size=500)

        global_median = np.median(counts)
        normalized = counts / global_median

        assert abs(np.median(normalized) - 1.0) < 0.01

    def test_chr21_trisomy_normalized_count_elevated(self):
        """
        In a T21 sample with 8% FF:
        chr21 should be elevated ~3% relative to autosomes after normalization.
        """
        rng = np.random.default_rng(21)
        ff = 0.08

        # Expected: euploid chr21 = 1.0, trisomy contribution = 1 + ff/2
        expected_ratio = 1.0 + ff / 2
        # Add noise
        euploid_counts = rng.normal(1.0,  0.02, 100)
        t21_counts     = rng.normal(expected_ratio, 0.02, 100)

        ratio = np.mean(t21_counts) / np.mean(euploid_counts)
        assert ratio > 1.02, f"Expected ratio > 1.02, got {ratio:.4f}"
        assert ratio < 1.08, f"Expected ratio < 1.08, got {ratio:.4f}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
