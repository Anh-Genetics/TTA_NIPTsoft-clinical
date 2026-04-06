# Troubleshooting Guide: TTA_NIPTsoft-clinical

**Document ID**: TTA-NIPTsoft-TRBL-001  
**Version**: 1.0.0  

---

## Quick Diagnosis

| Symptom | Likely Cause | Go To |
|---------|-------------|-------|
| Pipeline fails immediately | Samplesheet error | §1 |
| All samples get NO_CALL | Low read depth / bad QC | §2 |
| High no-call rate | Low fetal fraction batch | §3 |
| Alignment step fails | Missing/wrong reference index | §4 |
| GC correction fails | R package missing | §5 |
| Report generation fails | Missing upstream outputs | §6 |
| Unexpected HIGH RISK calls | Reference panel issue | §7 |
| Pipeline hangs | Resource limits | §8 |

---

## §1: Samplesheet Errors

### Error: `ERROR: Samplesheet validation failed`

**Diagnosis**: Run validation manually:
```bash
python3 bin/validate_samplesheet.py \
    --samplesheet your_samplesheet.csv \
    --schema schemas/samplesheet.schema.json \
    --out-csv /tmp/validated.csv \
    --out-report /tmp/report.json
cat /tmp/report.json | python3 -m json.tool
```

**Common causes and fixes**:

| Error | Fix |
|-------|-----|
| Duplicate sample_id | Ensure all sample_id values are unique |
| FASTQ file not found | Check absolute/relative paths; verify file exists |
| Invalid reported_sex | Use: `male`, `female`, or `unknown` |
| Empty samplesheet | Add at least one data row |
| Missing required column | Add `sample_id` and `fastq_path` columns |

### Error: `Path does not exist`
```
ERROR: FASTQ file not found: /path/to/sample.fastq.gz
```
- Check that the path is accessible from where Nextflow runs
- For cluster runs: ensure paths are on shared filesystem
- Use absolute paths

---

## §2: All Samples Getting NO_CALL

### Diagnosis

Check the QC decision files:
```bash
for f in results/qc/*.qc_decision.json; do
    echo "=== $f ===";
    python3 -c "
import json, sys
with open('$f') as fp:
    d = json.load(fp)
print('QC:', d['qc_status'])
for g in d.get('failed_gates', []):
    print(' FAILED:', g['code'], '-', g['message'],
          f'(observed={g[\"observed\"]}, threshold={g[\"threshold\"]})')
"
done
```

### Gate G01: Low Raw Reads

**Symptoms**: `G01_LOW_RAW_READS` in all samples

**Causes**:
- Sequencing run produced too few reads
- Wrong FASTQ file provided (e.g., undetermined reads)
- FASTQ file is empty or truncated

**Fix**:
```bash
# Check read counts
zcat sample.fastq.gz | wc -l  # divide by 4 = reads
# Check if file is valid
zcat sample.fastq.gz | head -8
```

### Gate G02: Low Unique Mapped Reads

**Symptoms**: `G02_LOW_UNIQUE_READS`

**Causes**:
- Wrong reference genome (mismatch with sequencing protocol)
- Contaminated library (non-human reads)
- Too-strict MAPQ threshold

**Fix**:
```bash
# Check alignment stats
cat results/qc/bam/SAMPLEID.bam_stats.txt
# Check top contaminating organisms (if samtools available)
samtools view -s 0.01 sample.filtered.bam | cut -f3 | sort | uniq -c | sort -rn | head
```

### Gate G04: Low Fetal Fraction

**Symptoms**: `G04_LOW_FF` — this is a clinical issue, not a technical error.

**Possible causes**:
- Very early gestational age (< 10 weeks)
- High maternal BMI (dilution effect)
- Sample storage/processing issues
- Maternal condition affecting cfDNA release

**Action**: 
1. Check gestational age in samplesheet
2. Review maternal BMI if available
3. Consider re-sampling at later gestational age
4. Consult clinical team

---

## §3: High No-Call Rate in a Batch

If > 20% of samples in a batch get no-calls:

1. **Check the sequencing run**: Review sequencer quality report
2. **Check batch controls**: Do internal controls pass?
3. **Check QC metrics across batch**:
   ```bash
   # Extract FF and reads for all samples in results
   python3 -c "
   import json, glob, os
   for f in sorted(glob.glob('results/qc/*.qc_decision.json')):
       with open(f) as fp:
           d = json.load(fp)
       m = d.get('qc_metrics', {})
       print(os.path.basename(f).replace('.qc_decision.json',''),
             'FF:', m.get('fetal_fraction'),
             'Map:', m.get('mapping_rate'),
             'Reads:', m.get('raw_reads'),
             'QC:', d['qc_status'])
   "
   ```
4. **Library QC**: Review fragment size distribution, adapter content

---

## §4: Alignment Step Fails

### Error: `Error running bowtie2`

**Check Bowtie2 index**:
```bash
ls -la assets/genome/hg19/
# Should see: hg19.1.bt2, hg19.2.bt2, hg19.3.bt2, hg19.4.bt2, hg19.rev.1.bt2, hg19.rev.2.bt2
bowtie2-inspect --summary assets/genome/hg19/hg19 2>&1 | head
```

**Check parameter**:
```bash
nextflow run main.nf --bowtie2_index /absolute/path/to/hg19
```

### Error: `[Error] index does not exist`
- Ensure Bowtie2 index base path is correct (without `.bt2` extension)
- Example: `--bowtie2_index /ref/hg19/hg19` (not `hg19.1.bt2`)

---

## §5: GC Correction Fails

### Error: `Error in gc_loess_correct.R`

**Check R dependencies**:
```R
# In R console
required <- c("optparse", "arrow", "jsonlite")
installed <- installed.packages()[,"Package"]
missing <- required[!required %in% installed]
if (length(missing)) cat("Missing packages:", paste(missing, collapse=", "), "\n")
```

**Install missing packages**:
```R
install.packages(c("optparse", "jsonlite"))
install.packages("arrow")
# or via BiocManager:
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
```

**Common errors**:

| Error | Cause | Fix |
|-------|-------|-----|
| `Error in loess()` | Too few valid bins | Check bin BED file; reduce min-gc/min-map thresholds |
| `Arrow package not found` | R arrow not installed | `install.packages("arrow")` |
| `Error reading BED file` | Format mismatch | Verify column order: chr, start, end, gc_frac |

---

## §6: Report Generation Fails

### Error: `Schema validation failed`

**Check report against schema**:
```bash
python3 -c "
import json, jsonschema
with open('results/reports/SAMPLEID_clinical_report.json') as f:
    report = json.load(f)
with open('schemas/report.schema.json') as f:
    schema = json.load(f)
jsonschema.validate(instance=report, schema=schema)
print('Schema validation PASSED')
"
```

### Missing upstream files

Report generation requires all upstream outputs. Check:
```bash
ls results/qc/SAMPLEID.qc_decision.json
ls results/calls/SAMPLEID.aneuploidy_calls.json
ls results/calls/SAMPLEID.sca_calls.json
```

---

## §7: Unexpected HIGH RISK Calls

### Systematic false positives (same chromosome across many samples)

This suggests a reference panel issue:
1. **Check reference panel distribution**:
   ```python
   import pandas as pd, pyarrow.parquet as pq
   ref = pq.read_table('assets/controls/euploid_reference_matrix.parquet').to_pandas()
   print(ref['cr_chr21'].describe())
   ```
2. If distribution is abnormal: re-build reference panel from fresh euploid runs
3. Check that reference panel was built with the same pipeline version

### Isolated unexpected call

1. Review QC metrics for the sample
2. Check for batch effects (compare Z-scores to other samples in run)
3. Review fetal fraction
4. Consider that it may be a true positive — recommend clinical follow-up

---

## §8: Pipeline Hangs or Runs Slowly

### Resource limits

```bash
# Check current resource allocation
cat nextflow.config | grep -A5 "process {"

# Increase resources
nextflow run main.nf \
    --max_cpus 32 \
    --max_memory 64.GB \
    -resume
```

### Nextflow work directory

```bash
# Check disk space
df -h work/
# Clean up old work directories (careful!)
nextflow clean -f -before last
```

### Check Nextflow log

```bash
cat .nextflow.log | grep ERROR
nextflow log last -f name,status,exit,realtime,cpus,memory,rss
```

---

## §9: Error Code Reference

| Code | Stage | Description |
|------|-------|-------------|
| E101 | Input validation | Samplesheet error |
| E102 | Input validation | FASTQ file not found |
| E201 | Alignment | Bowtie2 index error |
| E202 | Alignment | Samtools error |
| E301 | QC | Bin count empty |
| E401 | Statistics | GC correction failed |
| E501 | Calling | Reference panel error |
| E601 | Reporting | Schema validation failed |

---

## §10: Getting Help

1. Review the log: `cat .nextflow.log | tail -100`
2. Enable verbose logging: `nextflow run main.nf -with-trace -with-report`
3. Check issues: https://github.com/Anh-Genetics/TTA_NIPTsoft-clinical/issues
4. Include in bug reports:
   - Pipeline version (`cat nextflow.config | grep pipeline_version`)
   - Full error message
   - Sample QC metrics
   - Nextflow version (`nextflow -version`)
   - Operating system and Docker/Singularity version
