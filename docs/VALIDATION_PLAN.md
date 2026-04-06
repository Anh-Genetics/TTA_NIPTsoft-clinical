# Validation Plan: TTA_NIPTsoft-clinical

**Document ID**: TTA-NIPTsoft-VAL-001  
**Version**: 1.0.0  
**Status**: Draft — Required before clinical deployment  

---

## 1. Purpose

This document defines the validation strategy for TTA_NIPTsoft-clinical prior to clinical use. Each laboratory must complete this validation independently using their own patient cohort and sequencing platform.

---

## 2. Regulatory Context

This validation plan is designed to support:
- ISO 15189 accreditation requirements
- CLIA/CAP laboratory requirements (US)
- CE-IVD directive requirements (EU)
- Applicable national laboratory regulations

> This software is a **research tool** that requires laboratory-specific validation before clinical use.

---

## 3. Validation Studies Required

### 3.1 Study 1: Accuracy (Analytical Sensitivity/Specificity)

**Objective**: Demonstrate concordance between pipeline results and karyotype/FISH truth.

**Requirements**:
- Minimum cohort: ≥ 200 samples with known karyotype outcomes
  - T21: ≥ 30 positive cases
  - T18: ≥ 15 positive cases  
  - T13: ≥ 10 positive cases
  - Euploid controls: ≥ 100 samples
- Sequential samples processed as routine (not cherry-picked)
- Truth standard: diagnostic karyotype or FISH

**Acceptance Criteria (Tier 1)**:

| Target | Sensitivity | Specificity | PPV (at 1:1000 prevalence) |
|--------|------------|-------------|---------------------------|
| T21 | ≥ 99% | ≥ 99.5% | To be established |
| T18 | ≥ 97% | ≥ 99.5% | To be established |
| T13 | ≥ 90% | ≥ 99.5% | To be established |

**Note**: PPV depends on population prevalence and must be communicated in clinical reports.

### 3.2 Study 2: Reproducibility

**Objective**: Demonstrate consistent results across runs and instruments.

**Design**:
- **Intra-run**: Same sample processed 3× within one run
- **Inter-run**: Same sample processed across 5 independent runs
- **Inter-instrument**: If multiple sequencing instruments used

**Acceptance Criteria**:
- Z-score CV across replicates: ≤ 10%
- No discordant risk classifications between replicates

### 3.3 Study 3: Limit of Detection (LoD)

**Objective**: Define minimum fetal fraction for reliable detection.

**Design**:
- Use in vitro dilution series (known T21 samples diluted into euploid plasma)
- Test fetal fractions: 2%, 3%, 4%, 5%, 6%, 8%, 10%, 15%
- n ≥ 20 per fetal fraction level

**Acceptance Criteria**:
- LoD (95% detection): expected at FF ≥ 4% for T21
- Establish FF threshold for each target separately

### 3.4 Study 4: Interference Testing

**Objective**: Verify performance under challenging conditions.

**Test Conditions**:
- High maternal BMI (> 35 BMI) — tends to lower FF
- Multiple gestation pregnancies (twins)
- Prior blood transfusion (within 3 months)
- Maternal chromosomal abnormalities (if known)
- Low gestational age (< 11 weeks)

### 3.5 Study 5: No-Call Rate Characterization

**Objective**: Characterize no-call rate and reasons in routine population.

**Design**:
- Process ≥ 500 consecutive clinical samples
- Document all no-call reasons and frequency

**Acceptance Criteria**:
- Overall no-call rate: ≤ 5% (target ≤ 2%)
- Document distribution of no-call reasons

---

## 4. Extended Panel Validation (Tier 2)

### 4.1 SCA Validation

Additional requirements for sex chromosome aneuploidy (SCA) calls:
- Minimum: 20 positive cases each for 45,X and 47,XXY
- Minimum: 10 cases each for 47,XXX and 47,XYY
- Note: SCA PPV is generally lower than for T21/T18/T13

### 4.2 RAA Module Validation

The RAA scan module has **lower clinical performance** than Tier 1:
- Document positive predictive value in validation cohort
- Establish minimum evidence threshold for reporting
- Consider whether to report as "suspected" with clinical review required

### 4.3 Microdeletion Panel Validation

Microdeletion detection requires specific validation per syndrome region:
- Minimum: 10 positive cases per syndrome (where feasible)
- **Important**: PPV for microdeletions is typically low (< 50%) at common NIPT depths
- Laboratory must determine whether benefit outweighs harm from false positives

---

## 5. Reference Panel Validation

### 5.1 Euploid Reference Panel Requirements

| Parameter | Requirement |
|-----------|-------------|
| Minimum n | ≥ 100 euploid samples |
| Sex balance | Represent both male and female fetuses |
| GA range | Cover 10–22 weeks |
| Sequencing depth | Same protocol as routine testing |
| Processing | Same pipeline version as validation |

### 5.2 Reference Panel Monitoring

- Review reference panel statistics periodically (e.g., quarterly)
- Trigger re-validation if:
  - New sequencing reagent lot introduced
  - New instrument installed
  - Major pipeline update
  - QC metric drift detected

---

## 6. Bioinformatics Validation

### 6.1 Software Verification

Before first clinical use, verify:
- [ ] All unit tests pass: `pytest tests/unit/ -v`
- [ ] All integration tests pass: `pytest tests/integration/ -v`
- [ ] Pipeline runs to completion on positive control samples
- [ ] Pipeline runs to completion on euploid control samples
- [ ] No-call logic triggers correctly on low-FF samples
- [ ] Report schema validation passes
- [ ] Audit trail files generated correctly

### 6.2 Performance Benchmarking

Document on validation hardware:
- Run time per sample
- Memory requirements
- CPU requirements
- Storage requirements per run

### 6.3 Edge Case Testing

Test pipeline behavior with:
- Minimum read count (near QC threshold)
- Maximum read count
- Very high/low GC content samples
- All-no-call run
- Mixed pass/fail runs

---

## 7. Ongoing Quality Assurance

### 7.1 Internal Controls

Each run should include:
- ≥ 1 positive control (characterized aneuploidy sample or cell line)
- ≥ 1 negative control (euploid sample of known quality)
- Run fails QC if controls do not produce expected results

### 7.2 Performance Monitoring

Monitor monthly:
- No-call rate trend
- QC metric distributions (reads, mapping rate, FF)
- Z-score distributions for all targets
- Positive call rate (should be stable)

### 7.3 External Quality Assessment

Participate in external QA schemes when available:
- EMQN NIPT scheme (Europe)
- CAP proficiency testing (US)
- RCPA QAP programs (Australia)

---

## 8. Validation Documentation Requirements

Complete the following before clinical deployment:

- [ ] Validation study protocols approved by laboratory director
- [ ] Validation data collected and analyzed
- [ ] Validation report written and approved
- [ ] Acceptance criteria met or deviations documented
- [ ] QC thresholds locked in configuration
- [ ] SOP finalized and approved
- [ ] Staff training records completed
- [ ] Regulatory submission (if required)

---

## 9. Change Validation Requirements

When any of the following changes are made, partial or full re-validation is required:

| Change Type | Minimum Re-validation |
|-------------|----------------------|
| QC threshold change | Statistical justification + n ≥ 20 samples |
| Statistical model update | Accuracy study (full) |
| Reference genome change | Accuracy study (full) |
| New reference panel | Reproducibility + accuracy comparison |
| New library kit | Accuracy study (full) |
| New sequencer model | Reproducibility + accuracy comparison |

---

*This validation plan must be reviewed and adapted to local regulatory requirements by a qualified laboratory director.*
