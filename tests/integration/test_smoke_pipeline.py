"""
tests/integration/test_smoke_pipeline.py
Integration/smoke tests for core pipeline scripts.
Tests end-to-end data flow using mock data without requiring Nextflow.
"""

import sys
import os
import json
import subprocess
import tempfile
import shutil
import pytest
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
BIN_DIR   = REPO_ROOT / "bin"
TESTS_DIR = REPO_ROOT / "tests"


# ── Mock data factories ────────────────────────────────────────
def make_mock_norm_counts(tmp_dir: Path, sample_id: str,
                           t21_elevated: bool = False) -> Path:
    """Create a mock normalized counts Parquet file."""
    rng = np.random.default_rng(42)
    rows = []
    chrs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    for chr_name in chrs:
        n_bins = 40 if chr_name not in ("chrX", "chrY") else 20
        for i in range(n_bins):
            base = 1.0
            if t21_elevated and chr_name == "chr21":
                base = 1.05  # simulate trisomy 21 elevation (5% with 10% FF)

            rows.append({
                "sample_id":           sample_id,
                "chr":                 chr_name,
                "start":               i * 50000,
                "end":                 (i + 1) * 50000,
                "raw_count":           int(rng.normal(1000, 50)),
                "gc_corrected_count":  rng.normal(base, 0.002),
                "norm_count":          rng.normal(base, 0.002),
                "gc_frac":             rng.uniform(0.35, 0.65),
                "map_frac":            rng.uniform(0.7, 1.0),
                "is_masked":           False
            })

    df  = pd.DataFrame(rows)
    out = tmp_dir / f"{sample_id}.norm_counts.parquet"
    pq.write_table(pa.Table.from_pandas(df), str(out))
    return out


def make_mock_ff_metrics(tmp_dir: Path, sample_id: str,
                          ff: float = 0.08) -> Path:
    """Create a mock fetal fraction metrics JSON."""
    data = {
        "sample_id":    sample_id,
        "ff_seqff":     ff * 0.95,
        "ff_y":         ff * 1.05 if ff > 0 else None,
        "ff_final":     ff,
        "ff_confidence": 0.85,
        "ff_method":    "combined",
        "reported_sex": "unknown",
        "bam_metrics": {
            "total_reads":     12000000,
            "mapped_reads":    10800000,
            "duplicate_reads": 1200000,
            "mapping_rate":    0.90,
            "dup_rate":        0.11
        }
    }
    out = tmp_dir / f"{sample_id}.ff_metrics.json"
    out.write_text(json.dumps(data, indent=2))
    return out


def make_mock_euploid_reference(tmp_dir: Path) -> Path:
    """Create a mock euploid reference matrix Parquet.

    CR values are calibrated to match make_mock_norm_counts output:
    - Each autosome: 40 bins at norm_count ~1.0
    - NORM_CHRS (19 autosomes used as normalizers): 19 * 40 = 760 bins
    - Expected CR for any autosome ≈ 40 / 760 ≈ 0.0526
    """
    rng = np.random.default_rng(99)
    data = {}
    chrs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    # Calibrated base CRs matching the mock normalized counts structure
    # 40 bins per autosome, 19 normalizer autosomes → CR ≈ 40/760 = 0.05263
    AUTOSOME_BASE_CR = 40.0 / (19 * 40)  # ~0.0526
    for chr_name in chrs:
        if chr_name in ("chrX", "chrY"):
            base = 0.10  # sex chromosomes have different representation
        else:
            base = AUTOSOME_BASE_CR
        # Use tight reference std (~0.5% relative CV) as in real NIPT
        data[f"cr_{chr_name}"] = rng.normal(base, base * 0.005, 50)
    # Also add sex-chromosome specific
    data["cr_chrX_female"] = rng.normal(0.10, 0.003, 50)
    data["cr_chrX_male"]   = rng.normal(0.05, 0.002, 50)
    data["cr_chrY_male"]   = rng.normal(0.005, 0.0005, 50)
    data["cr_chrY_female"] = rng.normal(0.0001, 0.00005, 50)

    df  = pd.DataFrame(data)
    out = tmp_dir / "euploid_ref.parquet"
    pq.write_table(pa.Table.from_pandas(df), str(out))
    return out


def make_mock_fastp_json(tmp_dir: Path, sample_id: str) -> Path:
    """Create a mock fastp output JSON."""
    data = {
        "summary": {
            "before_filtering": {
                "total_reads": 12000000,
                "total_bases": 900000000,
                "q20_rate": 0.97,
                "q30_rate": 0.93,
                "gc_content": 0.52
            },
            "after_filtering": {
                "total_reads": 11800000,
                "total_bases": 885000000,
                "q20_rate": 0.99,
                "q30_rate": 0.96,
                "gc_content": 0.51
            }
        },
        "duplication": {"rate": 0.12}
    }
    out = tmp_dir / f"{sample_id}_fastp.json"
    out.write_text(json.dumps(data, indent=2))
    return out


def make_mock_bam_stats(tmp_dir: Path, sample_id: str) -> Path:
    """Create a mock samtools flagstat output."""
    content = (
        "12000000 + 0 in total (QC-passed reads + QC-failed reads)\n"
        "0 + 0 secondary\n"
        "0 + 0 supplementary\n"
        "1440000 + 0 duplicates\n"
        "10800000 + 0 mapped (90.00% : N/A)\n"
        "0 + 0 paired in sequencing\n"
    )
    out = tmp_dir / f"{sample_id}.bam_stats.txt"
    out.write_text(content)
    return out


def run_script(script: str, args: list) -> subprocess.CompletedProcess:
    """Run a bin/ Python script and return result."""
    cmd = [sys.executable, str(BIN_DIR / script)] + args
    return subprocess.run(cmd, capture_output=True, text=True)


# ── Tests ──────────────────────────────────────────────────────

class TestZscoreModelSmoke:
    """Smoke tests for zscore_model.py."""

    def test_euploid_sample_low_risk(self, tmp_path):
        sample_id = "SMOKE_EUPLOID"
        counts   = make_mock_norm_counts(tmp_path, sample_id, t21_elevated=False)
        ff       = make_mock_ff_metrics(tmp_path, sample_id, ff=0.08)
        ref      = make_mock_euploid_reference(tmp_path)
        out_json = tmp_path / f"{sample_id}.aneuploidy_calls.json"

        result = run_script("zscore_model.py", [
            "--norm-counts",        str(counts),
            "--ff-metrics",         str(ff),
            "--euploid-ref",        str(ref),
            "--sample-id",          sample_id,
            "--targets",            "chr21,chr18,chr13",
            "--zscore-high-risk",   "3.0",
            "--zscore-grey-zone",   "2.5",
            "--out",                str(out_json)
        ])
        assert result.returncode == 0, f"Script failed: {result.stderr}"
        assert out_json.exists()

        data = json.loads(out_json.read_text())
        assert data["sample_id"] == sample_id
        assert "chr21" in data["targets"]
        assert "chr18" in data["targets"]
        assert "chr13" in data["targets"]
        # Euploid sample should be LOW_RISK or borderline (not HIGH_RISK)
        assert data["targets"]["chr21"]["classification"] in ("LOW_RISK", "GREY_ZONE")

    def test_t21_sample_high_risk(self, tmp_path):
        """T21-elevated sample should be classified as HIGH_RISK."""
        sample_id = "SMOKE_T21"
        counts   = make_mock_norm_counts(tmp_path, sample_id, t21_elevated=True)
        ff       = make_mock_ff_metrics(tmp_path, sample_id, ff=0.10)
        ref      = make_mock_euploid_reference(tmp_path)
        out_json = tmp_path / f"{sample_id}.aneuploidy_calls.json"

        result = run_script("zscore_model.py", [
            "--norm-counts",        str(counts),
            "--ff-metrics",         str(ff),
            "--euploid-ref",        str(ref),
            "--sample-id",          sample_id,
            "--targets",            "chr21,chr18,chr13",
            "--zscore-high-risk",   "3.0",
            "--zscore-grey-zone",   "2.5",
            "--out",                str(out_json)
        ])
        assert result.returncode == 0, f"Script failed: {result.stderr}"
        data = json.loads(out_json.read_text())
        # With 4% elevation and tight reference, expect HIGH_RISK or at minimum GREY_ZONE
        assert data["targets"]["chr21"]["classification"] in ("HIGH_RISK", "GREY_ZONE"), \
            f"Expected elevated T21 classification, got: {data['targets']['chr21']}"

    def test_low_ff_triggers_no_call(self, tmp_path):
        """Sample with FF below threshold should get NO_CALL."""
        sample_id = "SMOKE_LOW_FF"
        counts   = make_mock_norm_counts(tmp_path, sample_id)
        ff       = make_mock_ff_metrics(tmp_path, sample_id, ff=0.02)
        ref      = make_mock_euploid_reference(tmp_path)
        out_json = tmp_path / f"{sample_id}.aneuploidy_calls.json"

        result = run_script("zscore_model.py", [
            "--norm-counts",        str(counts),
            "--ff-metrics",         str(ff),
            "--euploid-ref",        str(ref),
            "--sample-id",          sample_id,
            "--targets",            "chr21,chr18,chr13",
            "--zscore-high-risk",   "3.0",
            "--zscore-grey-zone",   "2.5",
            "--out",                str(out_json)
        ])
        assert result.returncode == 0
        data = json.loads(out_json.read_text())
        assert data["targets"]["chr21"]["classification"] == "NO_CALL"


class TestQCDecisionSmoke:
    """Smoke tests for qc_decision.py."""

    def test_passing_sample(self, tmp_path):
        sample_id    = "QC_PASS"
        fastp        = make_mock_fastp_json(tmp_path, sample_id)
        bam_stats    = make_mock_bam_stats(tmp_path, sample_id)
        ff_metrics   = make_mock_ff_metrics(tmp_path, sample_id, ff=0.09)

        # Create stub call files
        aneuploidy_calls = tmp_path / f"{sample_id}.aneuploidy_calls.json"
        aneuploidy_calls.write_text(json.dumps({"targets": {}, "ff_final": 0.09}))
        sca_calls = tmp_path / f"{sample_id}.sca_calls.json"
        sca_calls.write_text(json.dumps({"calls": {}, "ff_final": 0.09}))

        out_json = tmp_path / f"{sample_id}.qc_decision.json"

        result = run_script("qc_decision.py", [
            "--sample-id",         sample_id,
            "--fastp-json",        str(fastp),
            "--bam-stats",         str(bam_stats),
            "--ff-metrics",        str(ff_metrics),
            "--aneuploidy-calls",  str(aneuploidy_calls),
            "--sca-calls",         str(sca_calls),
            "--min-raw-reads",     "8000000",
            "--min-unique-reads",  "4000000",
            "--max-dup-rate",      "0.35",
            "--min-ff",            "0.04",
            "--min-mapping-rate",  "0.70",
            "--out",               str(out_json)
        ])
        assert result.returncode == 0, f"Script failed: {result.stderr}"
        data = json.loads(out_json.read_text())
        assert data["qc_status"] == "PASS"
        assert data["reportable"] is True
        assert data["failed_gates"] == []

    def test_failing_sample(self, tmp_path):
        """Sample with very low FF should fail QC."""
        sample_id    = "QC_FAIL"
        fastp        = make_mock_fastp_json(tmp_path, sample_id)
        bam_stats    = make_mock_bam_stats(tmp_path, sample_id)
        ff_metrics   = make_mock_ff_metrics(tmp_path, sample_id, ff=0.01)  # too low

        aneuploidy_calls = tmp_path / f"{sample_id}.aneuploidy_calls.json"
        aneuploidy_calls.write_text(json.dumps({"targets": {}, "ff_final": 0.01}))
        sca_calls = tmp_path / f"{sample_id}.sca_calls.json"
        sca_calls.write_text(json.dumps({"calls": {}, "ff_final": 0.01}))

        out_json = tmp_path / f"{sample_id}.qc_decision.json"

        result = run_script("qc_decision.py", [
            "--sample-id",         sample_id,
            "--fastp-json",        str(fastp),
            "--bam-stats",         str(bam_stats),
            "--ff-metrics",        str(ff_metrics),
            "--aneuploidy-calls",  str(aneuploidy_calls),
            "--sca-calls",         str(sca_calls),
            "--min-raw-reads",     "8000000",
            "--min-unique-reads",  "4000000",
            "--max-dup-rate",      "0.35",
            "--min-ff",            "0.04",
            "--min-mapping-rate",  "0.70",
            "--out",               str(out_json)
        ])
        assert result.returncode == 0
        data = json.loads(out_json.read_text())
        assert data["qc_status"] == "FAIL"
        assert data["reportable"] is False


class TestReportSmoke:
    """Smoke tests for make_report.py."""

    def test_report_generated(self, tmp_path):
        sample_id  = "RPT_SMOKE"
        ff_val     = 0.085

        # Build minimal input files
        qc_data = {
            "sample_id":       sample_id,
            "qc_status":       "PASS",
            "reportable":      True,
            "no_call_reasons": [],
            "failed_gates":    [],
            "qc_metrics": {
                "raw_reads":      12000000,
                "mapped_reads":   10800000,
                "dup_rate":       0.12,
                "mapping_rate":   0.90,
                "fetal_fraction": ff_val,
                "ff_method":      "combined"
            }
        }
        aneup_data = {
            "sample_id": sample_id,
            "model":     "robust_zscore_v1",
            "ff_final":  ff_val,
            "targets": {
                "chr21": {"chromosome": "chr21", "z_robust": 1.2, "classification": "LOW_RISK", "no_call_reason": None},
                "chr18": {"chromosome": "chr18", "z_robust": 0.8, "classification": "LOW_RISK", "no_call_reason": None},
                "chr13": {"chromosome": "chr13", "z_robust": -0.3, "classification": "LOW_RISK", "no_call_reason": None}
            }
        }
        sca_data = {
            "sample_id":    sample_id,
            "inferred_sex": "female",
            "reported_sex": "unknown",
            "ff_final":     ff_val,
            "calls": {
                "chrX": {"cr": 0.10, "z": 0.5, "classification": "LOW_RISK", "no_call_reason": None},
                "chrY": {"cr": 0.0001, "z": -1.2, "classification": "LOW_RISK", "no_call_reason": None}
            }
        }

        qc_file    = tmp_path / f"{sample_id}.qc_decision.json"
        aneup_file = tmp_path / f"{sample_id}.aneuploidy_calls.json"
        sca_file   = tmp_path / f"{sample_id}.sca_calls.json"
        schema_file = REPO_ROOT / "schemas" / "report.schema.json"

        qc_file.write_text(json.dumps(qc_data, indent=2))
        aneup_file.write_text(json.dumps(aneup_data, indent=2))
        sca_file.write_text(json.dumps(sca_data, indent=2))

        out_json  = tmp_path / f"{sample_id}_clinical_report.json"
        out_html  = tmp_path / f"{sample_id}_clinical_report.html"
        out_audit = tmp_path / f"{sample_id}.audit.json"

        result = run_script("make_report.py", [
            "--sample-id",        sample_id,
            "--run-id",           "RUN_20240101",
            "--qc-decision",      str(qc_file),
            "--aneuploidy-calls", str(aneup_file),
            "--sca-calls",        str(sca_file),
            "--pipeline-version", "1.0.0",
            "--genome-build",     "hg19",
            "--clinical-tier",    "tier1",
            "--config-version",   "1.0.0",
            "--report-schema",    str(schema_file),
            "--out-json",         str(out_json),
            "--out-html",         str(out_html),
            "--out-audit",        str(out_audit)
        ])
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        # Validate JSON structure
        report = json.loads(out_json.read_text())
        assert report["sample_id"]   == sample_id
        assert report["qc_summary"]["status"] == "PASS"
        assert "T21" in report["results"]
        assert "T18" in report["results"]
        assert "T13" in report["results"]
        assert report["results"]["T21"]["classification"] == "LOW_RISK"
        assert "disclaimer" in report

        # Validate HTML exists and has content
        html = out_html.read_text()
        assert len(html) > 1000
        assert sample_id in html
        assert "LOW RISK" in html

        # Validate audit file
        audit = json.loads(out_audit.read_text())
        assert audit["sample_id"]        == sample_id
        assert audit["pipeline_version"] == "1.0.0"


class TestValidateSamplesheetSmoke:
    """Smoke tests for validate_samplesheet.py."""

    def test_valid_samplesheet(self, tmp_path):
        # Create a dummy FASTQ file
        fastq = tmp_path / "sample1.fastq.gz"
        fastq.write_bytes(b"")  # empty file

        csv_content = f"sample_id,fastq_path,reported_sex,run_id\n"
        csv_content += f"SAMPLE001,{fastq},female,RUN001\n"

        ss_file  = tmp_path / "samplesheet.csv"
        out_csv  = tmp_path / "validated.csv"
        out_rpt  = tmp_path / "validation_report.json"
        ss_file.write_text(csv_content)

        schema = REPO_ROOT / "schemas" / "samplesheet.schema.json"

        result = run_script("validate_samplesheet.py", [
            "--samplesheet", str(ss_file),
            "--schema",      str(schema),
            "--out-csv",     str(out_csv),
            "--out-report",  str(out_rpt)
        ])
        assert result.returncode == 0, f"Validation failed: {result.stderr}"
        assert out_csv.exists()
        report = json.loads(out_rpt.read_text())
        assert report["status"] == "PASS"
        assert report["total_samples"] == 1

    def test_duplicate_sample_id_fails(self, tmp_path):
        fastq = tmp_path / "sample.fastq.gz"
        fastq.write_bytes(b"")

        csv_content = f"sample_id,fastq_path\n"
        csv_content += f"DUP001,{fastq}\n"
        csv_content += f"DUP001,{fastq}\n"  # duplicate

        ss_file = tmp_path / "samplesheet_dup.csv"
        ss_file.write_text(csv_content)

        result = run_script("validate_samplesheet.py", [
            "--samplesheet", str(ss_file),
            "--schema",      str(REPO_ROOT / "schemas" / "samplesheet.schema.json"),
            "--out-csv",     str(tmp_path / "out.csv"),
            "--out-report",  str(tmp_path / "report.json")
        ])
        assert result.returncode != 0  # Should fail

    def test_missing_required_field_fails(self, tmp_path):
        csv_content = "reported_sex\nfemale\n"  # missing sample_id and fastq_path
        ss_file = tmp_path / "incomplete.csv"
        ss_file.write_text(csv_content)

        result = run_script("validate_samplesheet.py", [
            "--samplesheet", str(ss_file),
            "--schema",      str(REPO_ROOT / "schemas" / "samplesheet.schema.json"),
            "--out-csv",     str(tmp_path / "out.csv"),
            "--out-report",  str(tmp_path / "report.json")
        ])
        assert result.returncode != 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
