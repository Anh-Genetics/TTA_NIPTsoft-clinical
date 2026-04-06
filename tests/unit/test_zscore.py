"""
tests/unit/test_zscore.py
Unit tests for zscore_model.py — aneuploidy calling statistics.
"""

import sys
import os
import json
import pytest
import numpy as np

# Add bin/ to path for direct imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../bin"))


# ── Import functions under test ──────────────────────────────
from zscore_model import (
    robust_zscore,
    ncv_zscore,
    classify_risk,
    compute_chromosome_representation,
    NORM_CHRS,
)


class TestRobustZscore:
    """Tests for robust_zscore() function."""

    def test_zero_for_mean_value(self):
        ref = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        assert robust_zscore(1.0, ref) == 0.0

    def test_positive_for_elevated_value(self):
        ref = np.array([1.0, 1.0, 1.0, 1.0, 1.1, 0.9, 1.02, 0.98, 1.0, 0.99])
        z = robust_zscore(1.5, ref)
        assert z > 3.0, f"Expected z > 3.0, got {z}"

    def test_negative_for_low_value(self):
        ref = np.array([1.0, 1.0, 1.0, 1.0, 1.1, 0.9, 1.02, 0.98, 1.0, 0.99])
        z = robust_zscore(0.5, ref)
        assert z < -3.0, f"Expected z < -3.0, got {z}"

    def test_returns_zero_when_mad_and_std_zero(self):
        ref = np.array([1.0, 1.0, 1.0])
        z = robust_zscore(1.0, ref)
        assert z == 0.0

    def test_handles_single_element_ref(self):
        ref = np.array([1.0])
        z = robust_zscore(2.0, ref)
        assert isinstance(z, float)

    def test_euploid_reference_gives_low_zscore(self):
        """Values within the reference distribution should give |Z| < 2."""
        rng = np.random.default_rng(42)
        ref = rng.normal(0.01, 0.0005, 100)
        value = np.median(ref)
        z = robust_zscore(value, ref)
        assert abs(z) < 1.0, f"Expected |z| < 1.0 for median value, got {z}"

    def test_t21_like_elevation(self):
        """Simulate T21-like elevation: ~3% increase in chr21 representation."""
        rng = np.random.default_rng(0)
        ref = rng.normal(0.01000, 0.00030, 200)
        trisomy_value = 0.01000 * 1.03 * (1 + 0.08)  # ~3% increase with 8% FF
        z = robust_zscore(trisomy_value, ref)
        assert z > 2.0, f"T21-like value should give Z > 2.0, got {z}"


class TestNcvZscore:
    """Tests for ncv_zscore() function."""

    def test_basic_positive_zscore(self):
        # Reference with non-zero std so Z-score is meaningful
        ref = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.1])
        z = ncv_zscore(2.0, ref)
        assert z > 0.0

    def test_symmetry(self):
        # Use symmetric reference so high/low deviations give similar |Z|
        rng = np.random.default_rng(5)
        ref = rng.normal(1.0, 0.05, 100)
        delta = 0.20  # same distance above and below mean
        z_high = ncv_zscore(np.mean(ref) + delta, ref)
        z_low  = ncv_zscore(np.mean(ref) - delta, ref)
        assert abs(abs(z_high) - abs(z_low)) < 0.5  # roughly symmetric

    def test_returns_zero_when_std_zero(self):
        ref = np.array([1.0, 1.0, 1.0])
        z = ncv_zscore(2.0, ref)
        assert z == 0.0


class TestClassifyRisk:
    """Tests for classify_risk() function."""

    def test_high_risk_above_threshold(self):
        assert classify_risk(3.5, high_thresh=3.0, grey_lo=2.5) == "HIGH_RISK"

    def test_exactly_at_high_threshold(self):
        assert classify_risk(3.0, high_thresh=3.0, grey_lo=2.5) == "HIGH_RISK"

    def test_grey_zone(self):
        assert classify_risk(2.7, high_thresh=3.0, grey_lo=2.5) == "GREY_ZONE"

    def test_low_risk_below_grey(self):
        assert classify_risk(1.0, high_thresh=3.0, grey_lo=2.5) == "LOW_RISK"

    def test_negative_zscore_is_low_risk(self):
        assert classify_risk(-1.0, high_thresh=3.0, grey_lo=2.5) == "LOW_RISK"

    def test_custom_thresholds(self):
        assert classify_risk(4.0, high_thresh=4.0, grey_lo=3.0) == "HIGH_RISK"
        assert classify_risk(3.5, high_thresh=4.0, grey_lo=3.0) == "GREY_ZONE"
        assert classify_risk(2.5, high_thresh=4.0, grey_lo=3.0) == "LOW_RISK"


class TestChrRepresentation:
    """Tests for compute_chromosome_representation()."""

    def _make_df(self):
        """Create a minimal mock normalized counts dataframe."""
        import pandas as pd
        data = []
        for chrom in ["chr21", "chr18", "chr13"] + NORM_CHRS[:5]:
            for i in range(10):
                data.append({
                    "chr": chrom, "start": i * 50000, "end": (i + 1) * 50000,
                    "norm_count": 1.0, "is_masked": False
                })
        return pd.DataFrame(data)

    def test_returns_float(self):
        df = self._make_df()
        cr = compute_chromosome_representation(df, "chr21", NORM_CHRS[:5])
        assert isinstance(cr, float)

    def test_euploid_representation_near_expected(self):
        """For equal bin counts, chr21 representation should be 10 bins / (5*10) bins = 0.2."""
        df = self._make_df()
        cr = compute_chromosome_representation(df, "chr21", NORM_CHRS[:5])
        assert abs(cr - 0.2) < 0.01, f"Expected ~0.2, got {cr}"

    def test_returns_nan_for_missing_chromosome(self):
        import pandas as pd
        df = self._make_df()
        cr = compute_chromosome_representation(df, "chrNONEXISTENT", NORM_CHRS[:5])
        assert cr == 0.0 or (cr != cr)  # 0.0 or NaN

    def test_elevated_chr21_increases_cr(self):
        import pandas as pd
        df = self._make_df()
        # Double the chr21 counts
        df.loc[df["chr"] == "chr21", "norm_count"] = 2.0
        cr_elevated = compute_chromosome_representation(df, "chr21", NORM_CHRS[:5])
        cr_normal   = compute_chromosome_representation(df, "chr18", NORM_CHRS[:5])
        assert cr_elevated > cr_normal

    def test_masked_bins_excluded(self):
        import pandas as pd
        df = self._make_df()
        # Mask all chr21 bins
        df.loc[df["chr"] == "chr21", "is_masked"] = True
        cr = compute_chromosome_representation(df, "chr21", NORM_CHRS[:5])
        assert cr == 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
