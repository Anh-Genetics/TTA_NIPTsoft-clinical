# Standard Operating Procedure: TTA_NIPTsoft-clinical

**Document ID**: TTA-NIPTsoft-SOP-001  
**Version**: 1.0.0  
**Effective Date**: 2024-01-01  
**Status**: Draft — Requires laboratory validation before clinical use  
**Author**: Anh-Genetics Bioinformatics Team

---

## 1. Purpose and Scope

This SOP describes the operational procedures for running the TTA_NIPTsoft-clinical pipeline for clinical NIPT analysis of Illumina single-end sequencing data.

**Scope**: Laboratory personnel responsible for NIPT data analysis.

---

## 2. Clinical Application

### 2.1 Intended Use

TTA_NIPTsoft-clinical analyzes cell-free DNA (cfDNA) from maternal blood to screen for:

**Tier 1 (mandatory):**
- Trisomy 21 (Down syndrome)
- Trisomy 18 (Edwards syndrome)
- Trisomy 13 (Patau syndrome)

**Tier 2 (extended, when validated):**
- Sex chromosome aneuploidies (Turner 45,X; Klinefelter 47,XXY; Triple X 47,XXX; 47,XYY)
- Rare autosomal aneuploidies (scan with strict flagging)
- Microdeletion panel (curated regions: 22q11.2, 1p36, 4p16.3, 5p15.2, 15q11-q13, others)

### 2.2 Limitations

- NIPT is a **screening test only** — not a diagnostic test
- Results cannot be used to diagnose chromosomal conditions
- High-risk results **must** be confirmed by invasive diagnostic testing
- Performance may be affected by: low fetal fraction, high maternal BMI, multiple gestation, maternal chromosomal abnormalities, recent blood transfusion

---

## 3. Pre-Analytical Requirements

### 3.1 Sample Requirements

| Parameter | Requirement |
|-----------|-------------|
| Sample type | Maternal plasma cfDNA |
| Gestational age | ≥ 10 weeks recommended |
| Library type | Illumina single-end, 35–75 bp |
| Sequencing depth | ≥ 10M reads/sample (tier1), ≥ 12M (tier2) |

### 3.2 Samplesheet Preparation

1. Create a CSV samplesheet with required columns: `sample_id`, `fastq_path`
2. Optional columns: `reported_sex`, `run_id`, `batch`, `ga_weeks`
3. Validate the samplesheet: `python3 bin/validate_samplesheet.py --samplesheet <file> --schema schemas/samplesheet.schema.json`
4. Ensure `sample_id` values are unique
5. Verify all FASTQ paths exist and are readable

### 3.3 Reference Setup

1. Place reference files in `assets/genome/`:
   - `hg19/` — Bowtie2 index files
   - `hg19_50kb_bins.bed` — Genomic bin definitions
   - `hg19_50kb_gc.bed` — GC content per bin
   - `hg19_50kb_mappability.bed` — Mappability per bin
   - `hg19_blacklist.bed` — ENCODE blacklist regions
2. Place control files in `assets/controls/`:
   - `euploid_reference_matrix.parquet` — Validated euploid reference panel

---

## 4. Analysis Procedure

### 4.1 Standard Run (Tier 1)

```bash
nextflow run main.nf \
    --samplesheet /path/to/samplesheet.csv \
    --clinical_tier tier1 \
    -profile docker,illumina_se,clinical_tier1 \
    --outdir /path/to/results \
    -resume
```

### 4.2 Extended Run (Tier 2)

```bash
nextflow run main.nf \
    --samplesheet /path/to/samplesheet.csv \
    --clinical_tier tier2 \
    --enable_raa true \
    --enable_microdeletion true \
    -profile docker,illumina_se,clinical_tier2 \
    --outdir /path/to/results \
    -resume
```

### 4.3 Run Monitoring

1. Monitor execution: `nextflow log -f name,status,exit,duration`
2. Check for errors: `grep ERROR .nextflow.log`
3. Review pipeline report: `results/pipeline_info/execution_report.html`

---

## 5. QC Evaluation

### 5.1 Run-level QC

Before releasing any results, verify:
1. All samples completed all pipeline steps
2. Internal controls (if used) pass expected ranges
3. No batch-level anomalies in QC metrics

### 5.2 Sample-level QC Gates

Review `results/qc/SAMPLEID.qc_decision.json` for each sample:

| Gate Code | Description | Action if Failed |
|-----------|-------------|-----------------|
| G01_LOW_RAW_READS | Raw reads below minimum | Repeat sequencing |
| G02_LOW_UNIQUE_READS | Mapped reads below minimum | Check alignment, repeat |
| G03_HIGH_DUP_RATE | High duplication rate | Review library prep |
| G04_LOW_FF | Fetal fraction < 4% | Review gestational age, resample |
| G05_LOW_MAPPING_RATE | Poor alignment | Check sample/reference integrity |
| G06_MISSING_FF | FF could not be estimated | Review data quality |

### 5.3 No-Call Criteria

A no-call result is issued when:
- Any mandatory QC gate fails
- Statistical confidence is insufficient
- Data quality indicators are outside validated range

**No-call results require clinical review and patient management decision.**

---

## 6. Result Interpretation

### 6.1 Risk Classifications

| Classification | Meaning | Clinical Action |
|----------------|---------|-----------------|
| HIGH RISK | Elevated chromosome representation consistent with aneuploidy | Genetic counseling; offer diagnostic testing |
| LOW RISK | Normal chromosome representation | Communicate result with screening test limitations |
| BORDERLINE | Intermediate Z-score (grey zone) | Clinical review; consider repeat testing |
| NO CALL | Insufficient data or QC failure | Document reason; consider repeat sample |
| SUSPECTED RAA (tier2) | Potential rare autosomal aneuploidy | Clinical review; genetic counseling |

### 6.2 Result Review Workflow

1. Automated report generated: `results/reports/SAMPLEID_clinical_report.json`
2. Bioinformatician reviews QC metrics
3. Clinical review by qualified personnel for all HIGH RISK and BORDERLINE results
4. Final report authorized and released to requesting clinician
5. All NO CALL results reviewed and documented with reason

---

## 7. Documentation and Records

### 7.1 Required Records

For each run, retain:
1. Samplesheet (original)
2. Pipeline execution log (`.nextflow.log`)
3. Pipeline report (`execution_report.html`)
4. Per-sample QC decisions (`*.qc_decision.json`)
5. Per-sample clinical reports (`*_clinical_report.json`)
6. Audit trail files (`*.audit.json`)
7. Software versions file

### 7.2 Audit Trail

Each report includes an audit file (`SAMPLEID.audit.json`) containing:
- Pipeline version and configuration version
- Input file checksums (MD5)
- Timestamp
- Environment information

---

## 8. Configuration Change Control

> **IMPORTANT**: Clinical QC thresholds and calling model thresholds are validated operating points. Any changes require formal change control review.

Controlled parameters (in `conf/clinical_tier1.config`, `conf/clinical_tier2.config`):
- `min_raw_reads`, `min_unique_reads`, `max_dup_rate`
- `min_fetal_fraction`, `min_mapping_rate`
- `zscore_high_risk`, `zscore_grey_zone_lo`
- `sca_zscore_threshold`, `raa_zscore_threshold`, `microdeletion_zscore`

Changes to controlled parameters require:
1. Analytical validation of the new threshold
2. Review and approval by laboratory director
3. Documentation in change control log
4. Update of `config_version` parameter
5. Communication to clinical staff

---

## 9. Troubleshooting

See `docs/TROUBLESHOOTING.md` for detailed troubleshooting guidance.

---

## 10. Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0.0 | 2024-01-01 | Anh-Genetics | Initial draft |
