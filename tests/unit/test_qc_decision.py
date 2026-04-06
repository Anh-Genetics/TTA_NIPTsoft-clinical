"""
tests/unit/test_qc_decision.py
Unit tests for qc_decision.py — QC gate logic.
"""

import sys
import os
import json
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../bin"))

from qc_decision import (
    evaluate_qc_gates,
    parse_fastp_json,
    parse_bam_stats,
    QC_GATE_CODES,
)


# ── Fixtures ─────────────────────────────────────────────────
def make_passing_fastp():
    return {
        "raw_reads":        12000000,
        "clean_reads":      11500000,
        "q20_rate":         0.98,
        "q30_rate":         0.95,
        "gc_content":       0.52,
        "duplication_rate": 0.15
    }

def make_passing_bam():
    return {
        "total_reads":     11500000,
        "mapped_reads":    10000000,
        "duplicate_reads": 1500000,
        "mapping_rate":    0.87,
        "dup_rate":        0.15
    }

def make_passing_ff():
    return {
        "ff_final":    0.08,
        "ff_method":   "combined",
        "ff_seqff":    0.075,
        "ff_y":        0.085
    }

THRESHOLDS = dict(
    min_raw_reads    = 8000000,
    min_unique_reads = 4000000,
    max_dup_rate     = 0.35,
    min_ff           = 0.04,
    min_mapping_rate = 0.70
)


class TestEvaluateQCGates:
    """Tests for evaluate_qc_gates()."""

    def test_all_passing_returns_pass(self):
        status, gates, metrics = evaluate_qc_gates(
            make_passing_fastp(), make_passing_bam(), make_passing_ff(),
            **THRESHOLDS
        )
        assert status == "PASS"
        assert gates == []

    def test_low_raw_reads_fails_g01(self):
        fastp = make_passing_fastp()
        fastp["raw_reads"] = 1000
        bam = make_passing_bam()
        bam["total_reads"] = 1000
        status, gates, _ = evaluate_qc_gates(
            fastp, bam, make_passing_ff(), **THRESHOLDS
        )
        assert status == "FAIL"
        codes = [g["code"] for g in gates]
        assert "G01_LOW_RAW_READS" in codes

    def test_low_unique_reads_fails_g02(self):
        bam = make_passing_bam()
        bam["mapped_reads"] = 100000
        bam["mapping_rate"] = 0.01
        status, gates, _ = evaluate_qc_gates(
            make_passing_fastp(), bam, make_passing_ff(), **THRESHOLDS
        )
        assert status == "FAIL"
        codes = [g["code"] for g in gates]
        assert "G02_LOW_UNIQUE_READS" in codes

    def test_high_dup_rate_fails_g03(self):
        bam = make_passing_bam()
        bam["dup_rate"] = 0.60
        status, gates, _ = evaluate_qc_gates(
            make_passing_fastp(), bam, make_passing_ff(), **THRESHOLDS
        )
        assert status == "FAIL"
        codes = [g["code"] for g in gates]
        assert "G03_HIGH_DUP_RATE" in codes

    def test_low_ff_fails_g04(self):
        ff = make_passing_ff()
        ff["ff_final"] = 0.02
        status, gates, _ = evaluate_qc_gates(
            make_passing_fastp(), make_passing_bam(), ff, **THRESHOLDS
        )
        assert status == "FAIL"
        codes = [g["code"] for g in gates]
        assert "G04_LOW_FF" in codes

    def test_missing_ff_fails_g06(self):
        ff = make_passing_ff()
        ff["ff_final"] = None
        status, gates, _ = evaluate_qc_gates(
            make_passing_fastp(), make_passing_bam(), ff, **THRESHOLDS
        )
        assert status == "FAIL"
        codes = [g["code"] for g in gates]
        assert "G06_MISSING_FF" in codes

    def test_low_mapping_rate_fails_g05(self):
        bam = make_passing_bam()
        bam["mapping_rate"] = 0.30
        status, gates, _ = evaluate_qc_gates(
            make_passing_fastp(), bam, make_passing_ff(), **THRESHOLDS
        )
        assert status == "FAIL"
        codes = [g["code"] for g in gates]
        assert "G05_LOW_MAPPING_RATE" in codes

    def test_multiple_failures_reported(self):
        """Both low reads and low FF should both be flagged."""
        fastp = make_passing_fastp()
        fastp["raw_reads"] = 1000
        bam = make_passing_bam()
        bam["total_reads"] = 1000
        bam["mapped_reads"] = 500
        ff = {"ff_final": 0.01, "ff_method": "seqff"}
        status, gates, _ = evaluate_qc_gates(
            fastp, bam, ff, **THRESHOLDS
        )
        assert status == "FAIL"
        assert len(gates) >= 2

    def test_exact_threshold_boundary_raw_reads(self):
        """Exactly at minimum threshold should PASS."""
        fastp = make_passing_fastp()
        fastp["raw_reads"] = 8000000
        bam = make_passing_bam()
        bam["total_reads"] = 8000000
        status, gates, _ = evaluate_qc_gates(
            fastp, bam, make_passing_ff(), **THRESHOLDS
        )
        # Exactly at limit: should NOT fail G01
        codes = [g["code"] for g in gates]
        assert "G01_LOW_RAW_READS" not in codes

    def test_metrics_returned(self):
        """Verify all expected metrics are returned."""
        _, _, metrics = evaluate_qc_gates(
            make_passing_fastp(), make_passing_bam(), make_passing_ff(),
            **THRESHOLDS
        )
        assert "fetal_fraction" in metrics
        assert "mapping_rate" in metrics
        assert "raw_reads" in metrics
        assert "mapped_reads" in metrics
        assert "dup_rate" in metrics


class TestParseFastpJson:
    """Tests for parse_fastp_json() using temp files."""

    def test_parses_correctly(self, tmp_path):
        fastp_output = {
            "summary": {
                "before_filtering": {"total_reads": 10000000},
                "after_filtering":  {
                    "total_reads": 9800000,
                    "q20_rate":    0.98,
                    "q30_rate":    0.95,
                    "gc_content":  0.51
                }
            },
            "duplication": {"rate": 0.12}
        }
        f = tmp_path / "fastp.json"
        f.write_text(json.dumps(fastp_output))
        result = parse_fastp_json(str(f))
        assert result["raw_reads"]        == 10000000
        assert result["clean_reads"]      == 9800000
        assert result["q20_rate"]         == 0.98
        assert result["duplication_rate"] == 0.12

    def test_returns_defaults_for_missing_file(self):
        result = parse_fastp_json("/nonexistent/file.json")
        assert result["raw_reads"] == 0
        assert result["q20_rate"]  == 0.0


class TestParseBamStats:
    """Tests for parse_bam_stats()."""

    def test_parses_flagstat(self, tmp_path):
        flagstat = """10000000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
1500000 + 0 duplicates
8700000 + 0 mapped (87.00% : N/A)
0 + 0 paired in sequencing
"""
        f = tmp_path / "bam.stats.txt"
        f.write_text(flagstat)
        result = parse_bam_stats(str(f))
        assert result["total_reads"]     == 10000000
        assert result["mapped_reads"]    == 8700000
        assert result["duplicate_reads"] == 1500000
        assert abs(result["mapping_rate"] - 0.87) < 0.001

    def test_returns_defaults_for_missing_file(self):
        result = parse_bam_stats("/nonexistent/bam.stats.txt")
        assert result["total_reads"]  == 0
        assert result["mapping_rate"] == 0.0


class TestQCGateCodes:
    """Ensure all gate codes have descriptions."""

    def test_all_codes_have_descriptions(self):
        expected_codes = [
            "G01_LOW_RAW_READS", "G02_LOW_UNIQUE_READS", "G03_HIGH_DUP_RATE",
            "G04_LOW_FF", "G05_LOW_MAPPING_RATE", "G06_MISSING_FF"
        ]
        for code in expected_codes:
            assert code in QC_GATE_CODES, f"Missing description for code: {code}"
            assert len(QC_GATE_CODES[code]) > 10


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
