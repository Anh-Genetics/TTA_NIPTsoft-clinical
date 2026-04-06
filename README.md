# TTA_NIPTsoft-clinical

**Clinical-grade Non-Invasive Prenatal Testing (NIPT) analysis pipeline for Illumina single-end sequencing data.**

> ⚠️ **CLINICAL DISCLAIMER**: NIPT is a **screening test and is not diagnostic**. A positive (high-risk) result does not confirm a chromosomal abnormality, and a negative (low-risk) result does not exclude one. All high-risk results should be confirmed by diagnostic testing (e.g., amniocentesis or chorionic villus sampling) following genetic counseling. This software must be validated by each laboratory prior to clinical use in accordance with applicable regulations.

---

## Overview

TTA_NIPTsoft-clinical is a production-ready Nextflow DSL2 pipeline for clinical NIPT analysis supporting:

| Target                     | Tier   | Classification               |
|---------------------------|--------|------------------------------|
| Trisomy 21 (T21)          | Tier 1 | High risk / Low risk / No-call |
| Trisomy 18 (T18)          | Tier 1 | High risk / Low risk / No-call |
| Trisomy 13 (T13)          | Tier 1 | High risk / Low risk / No-call |
| Sex chromosome aneuploidies (SCAs) | Tier 2 | High risk / Low risk / No-call |
| Rare autosomal aneuploidies (RAAs) | Tier 2 | Suspected RAA / Low risk / No-call |
| Microdeletion panel        | Tier 2 | High risk / Low risk / No-call |

## Architecture

```
FASTQ (SE)
   │
   ├─► VALIDATE_SAMPLESHEET     — Schema validation of input CSV
   │
   ├─► FASTP_QC                 — Read QC, trimming, adapter removal
   │
   ├─► BOWTIE2_ALIGN            — Alignment to reference genome (hg19)
   │
   ├─► BAM_FILTER               — MAPQ filter, blacklist removal, dup marking
   │
   ├─► BIN_COUNT                — Read counting in 50kb genomic bins
   │
   ├─► GC_CORRECT               — LOESS GC/mappability bias correction
   │
   ├─► FF_ESTIMATE              — Fetal fraction (SeqFF + Y-based)
   │
   ├─► ANEUPLOIDY_CALL          — T21/T18/T13 robust Z-score calling
   │
   ├─► SCA_CALL                 — chrX/chrY SCA calling
   │
   ├─► RAA_SCAN (tier2)         — Autosome-wide RAA scan
   │
   ├─► MICRODELETION_CALL (tier2) — Panel-constrained microdeletion calling
   │
   ├─► QC_GATE                  — Unified QC gate + no-call engine
   │
   └─► MAKE_REPORT              — JSON/HTML clinical report + audit trail
```

## Requirements

- **Nextflow** ≥ 23.04.0
- **Docker** or **Singularity** (recommended for production)
- **Python** ≥ 3.10 with: numpy, pandas, scipy, pyarrow, jsonschema
- **R** ≥ 4.3 with: optparse, arrow, jsonlite
- Reference files: hg19 genome, Bowtie2 index, GC/mappability BEDs, euploid reference matrix

## Quick Start

### 1. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### 2. Prepare samplesheet (CSV)
```csv
sample_id,fastq_path,reported_sex,run_id,batch,ga_weeks
SAMPLE001,/data/fastq/S001.fastq.gz,unknown,RUN20240101,batch_A,12
SAMPLE002,/data/fastq/S002.fastq.gz,female,RUN20240101,batch_A,14
```

### 3. Run Tier 1 (T21/T18/T13)
```bash
nextflow run main.nf \
    --samplesheet samplesheet.csv \
    --clinical_tier tier1 \
    --bowtie2_index /ref/hg19/hg19 \
    --euploid_reference /ref/controls/euploid_ref.parquet \
    -profile docker,illumina_se,clinical_tier1
```

### 4. Run Tier 2 (full extended panel)
```bash
nextflow run main.nf \
    --samplesheet samplesheet.csv \
    --clinical_tier tier2 \
    --enable_raa true \
    --enable_microdeletion true \
    -profile docker,illumina_se,clinical_tier2
```

## Configuration Profiles

| Profile | Description |
|---------|-------------|
| `illumina_se` | Illumina single-end read settings |
| `clinical_tier1` | T21/T18/T13 calling thresholds |
| `clinical_tier2` | Extended panel (SCA/RAA/microdeletion) |
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers |
| `test` | Minimal test run with mock data |

## Input Samplesheet Schema

Defined in `schemas/samplesheet.schema.json`.

**Required fields:**
- `sample_id` — unique identifier (alphanumeric, dash, underscore)
- `fastq_path` — path to input FASTQ.gz file

**Optional fields:**
- `reported_sex` — `male` / `female` / `unknown`
- `run_id`, `batch`, `ga_weeks`, `maternal_age`

## Outputs

```
results/
├── reports/
│   ├── SAMPLEID_clinical_report.json   ← Machine-readable report
│   └── SAMPLEID_clinical_report.html   ← Human-readable report
├── calls/
│   ├── SAMPLEID.aneuploidy_calls.json
│   ├── SAMPLEID.sca_calls.json
│   └── SAMPLEID.raa_calls.json (tier2)
├── qc/
│   ├── SAMPLEID.qc_decision.json
│   ├── fastp/ (per-sample QC)
│   └── bam/   (alignment stats)
├── counts/
│   └── SAMPLEID.norm_counts.parquet
├── stats/
│   └── SAMPLEID.ff_metrics.json
└── audit/
    ├── SAMPLEID.audit.json
    └── validation_report.json
```

## Clinical Report Structure

```json
{
  "schema_version": "1.0",
  "pipeline_version": "1.0.0",
  "genome_build": "hg19",
  "clinical_tier": "tier1",
  "sample_id": "SAMPLE001",
  "run_id": "RUN20240101",
  "qc_summary": {
    "status": "PASS",
    "reportable": true,
    "fetal_fraction": 0.082,
    "mapping_rate": 0.895,
    "raw_reads": 12500000
  },
  "results": {
    "T21": { "classification": "LOW_RISK", "clinical_risk": "LOW RISK", "z_score": 1.2 },
    "T18": { "classification": "LOW_RISK", "clinical_risk": "LOW RISK", "z_score": 0.8 },
    "T13": { "classification": "LOW_RISK", "clinical_risk": "LOW RISK", "z_score": -0.3 }
  },
  "inferred_sex": "female",
  "disclaimer": "NIPT is a screening test and is not diagnostic..."
}
```

## QC Gates

A result is reportable only when all QC gates pass:

| Gate | Metric | Default Threshold |
|------|--------|-------------------|
| G01 | Raw reads | ≥ 10M (tier1) |
| G02 | Unique mapped reads | ≥ 5M |
| G03 | Duplication rate | ≤ 30% |
| G04 | Fetal fraction | ≥ 4% |
| G05 | Mapping rate | ≥ 75% |

Thresholds are validated operating points and must not be changed without a formal change control review. See `docs/SOP.md`.

## Statistical Model

### Aneuploidy Calling (T21/T18/T13)
1. Compute **chromosome representation** (CR) = sum(target chr bins) / sum(normalizer chr bins)
2. Compare CR to euploid reference distribution using **robust Z-score** (MAD-based)
3. Apply decision thresholds:
   - Z ≥ 3.0 → **HIGH RISK**
   - 2.5 ≤ Z < 3.0 → **BORDERLINE** (recommend repeat)
   - Z < 2.5 → **LOW RISK**

### Fetal Fraction Estimation
- **SeqFF-inspired**: uses autosomal bin profile variance (sex-independent)
- **Y-based**: chrY representation × 2 (male fetuses only)
- **Combined**: weighted average when both methods available

### SCA Calling
- Separate sex-specific reference distributions for chrX/chrY
- Infers fetal sex from chrY representation
- Separate Z-score model for sex chromosome analysis

## Testing

```bash
# Install test dependencies
pip install pytest numpy pandas scipy pyarrow jsonschema

# Unit tests
pytest tests/unit/ -v

# Integration/smoke tests
python3 tests/data/generate_mock_data.py
pytest tests/integration/ -v

# All tests
pytest tests/ -v --tb=short
```

## Development

```bash
git clone https://github.com/Anh-Genetics/TTA_NIPTsoft-clinical
cd TTA_NIPTsoft-clinical

# Run with test profile (uses mock data)
nextflow run main.nf -profile test,docker
```

## Clinical Interpretation Boundaries

1. **NIPT is screening, not diagnostic** — all high-risk results require confirmatory testing
2. **Tier 2 modules** (RAA, microdeletion) have lower analytical performance than Tier 1
3. **No-call results** require clinical review; causes include low fetal fraction, insufficient reads, or technical failure
4. **Fetal fraction < 4%** triggers no-call for all targets
5. **Borderline Z-scores** (2.5–3.0) should prompt clinical discussion

## Version History

| Version | Date | Notes |
|---------|------|-------|
| 1.0.0 | 2024 | Initial release — T21/T18/T13 + SCA + RAA + microdeletion |

## License

This software is intended for research use and requires laboratory-specific validation prior to clinical deployment. See `docs/VALIDATION_PLAN.md`.

## Citation

If you use this pipeline in research, please cite:
> TTA_NIPTsoft-clinical v1.0.0. Anh-Genetics. https://github.com/Anh-Genetics/TTA_NIPTsoft-clinical