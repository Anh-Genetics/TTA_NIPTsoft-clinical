#!/usr/bin/env python3
"""
qc_decision.py
Unified QC gate and no-call decision engine.

Evaluates all QC metrics and determines:
1. Overall QC status (PASS/FAIL)
2. No-call reasons (if any)
3. Whether each target result should be reported or suppressed

Decision logic (all gates must pass for a reportable result):
- G01: Minimum raw reads
- G02: Minimum unique mapped reads
- G03: Maximum duplication rate
- G04: Minimum fetal fraction
- G05: Minimum mapping rate
"""

import argparse
import json
import sys
from datetime import datetime, timezone


QC_GATE_CODES = {
    "G01_LOW_RAW_READS":     "Raw read count below minimum threshold",
    "G02_LOW_UNIQUE_READS":  "Unique mapped reads below minimum threshold",
    "G03_HIGH_DUP_RATE":     "Duplication rate exceeds maximum threshold",
    "G04_LOW_FF":            "Fetal fraction below minimum threshold",
    "G05_LOW_MAPPING_RATE":  "Mapping rate below minimum threshold",
    "G06_MISSING_FF":        "Fetal fraction could not be estimated",
    "G07_MISSING_BAM_STATS": "BAM statistics could not be parsed"
}


def parse_fastp_json(path: str) -> dict:
    """Extract summary metrics from fastp JSON output."""
    try:
        with open(path) as f:
            data = json.load(f)
        summary = data.get("summary", {})
        before  = summary.get("before_filtering", {})
        after   = summary.get("after_filtering",  {})
        return {
            "raw_reads":   before.get("total_reads", 0),
            "clean_reads": after.get("total_reads", 0),
            "q20_rate":    after.get("q20_rate",  0.0),
            "q30_rate":    after.get("q30_rate",  0.0),
            "gc_content":  after.get("gc_content", 0.0),
            "duplication_rate": data.get("duplication", {}).get("rate", 0.0)
        }
    except (FileNotFoundError, json.JSONDecodeError, KeyError):
        return {
            "raw_reads": 0, "clean_reads": 0,
            "q20_rate": 0.0, "q30_rate": 0.0,
            "gc_content": 0.0, "duplication_rate": 0.0
        }


def parse_bam_stats(path: str) -> dict:
    """Parse samtools flagstat output."""
    metrics = {
        "total_reads":     0,
        "mapped_reads":    0,
        "duplicate_reads": 0,
        "mapping_rate":    0.0,
        "dup_rate":        0.0
    }
    try:
        with open(path) as f:
            for line in f:
                if "in total" in line:
                    metrics["total_reads"]  = int(line.split()[0])
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


def parse_ff_metrics(path: str) -> dict:
    """Parse fetal fraction metrics JSON."""
    try:
        with open(path) as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def evaluate_qc_gates(
    fastp: dict, bam: dict, ff: dict,
    min_raw_reads: int, min_unique_reads: int,
    max_dup_rate: float, min_ff: float, min_mapping_rate: float
) -> tuple[str, list, dict]:
    """
    Evaluate all QC gates.
    Returns: (status, failed_gates, all_metrics)
    """
    failed_gates = []
    metrics = {}

    # Extract values
    raw_reads    = fastp.get("raw_reads", bam.get("total_reads", 0))
    mapped_reads = bam.get("mapped_reads", 0)
    dup_rate     = bam.get("dup_rate", fastp.get("duplication_rate", 0.0))
    map_rate     = bam.get("mapping_rate", 0.0)
    ff_val       = ff.get("ff_final")

    metrics = {
        "raw_reads":      raw_reads,
        "mapped_reads":   mapped_reads,
        "dup_rate":       round(dup_rate, 4),
        "mapping_rate":   round(map_rate, 4),
        "fetal_fraction": round(ff_val, 4) if ff_val is not None else None,
        "ff_method":      ff.get("ff_method", "unknown"),
        "q20_rate":       round(fastp.get("q20_rate", 0.0), 4),
        "q30_rate":       round(fastp.get("q30_rate", 0.0), 4)
    }

    # Gate G01: Raw reads
    if raw_reads < min_raw_reads:
        failed_gates.append({
            "code": "G01_LOW_RAW_READS",
            "message": QC_GATE_CODES["G01_LOW_RAW_READS"],
            "observed": raw_reads,
            "threshold": min_raw_reads
        })

    # Gate G02: Unique mapped reads
    if mapped_reads < min_unique_reads:
        failed_gates.append({
            "code": "G02_LOW_UNIQUE_READS",
            "message": QC_GATE_CODES["G02_LOW_UNIQUE_READS"],
            "observed": mapped_reads,
            "threshold": min_unique_reads
        })

    # Gate G03: Duplication rate
    if dup_rate > max_dup_rate:
        failed_gates.append({
            "code": "G03_HIGH_DUP_RATE",
            "message": QC_GATE_CODES["G03_HIGH_DUP_RATE"],
            "observed": round(dup_rate, 4),
            "threshold": max_dup_rate
        })

    # Gate G04: Fetal fraction
    if ff_val is None:
        failed_gates.append({
            "code": "G06_MISSING_FF",
            "message": QC_GATE_CODES["G06_MISSING_FF"],
            "observed": None,
            "threshold": min_ff
        })
    elif ff_val < min_ff:
        failed_gates.append({
            "code": "G04_LOW_FF",
            "message": QC_GATE_CODES["G04_LOW_FF"],
            "observed": round(ff_val, 4),
            "threshold": min_ff
        })

    # Gate G05: Mapping rate
    if map_rate < min_mapping_rate and mapped_reads > 0:
        failed_gates.append({
            "code": "G05_LOW_MAPPING_RATE",
            "message": QC_GATE_CODES["G05_LOW_MAPPING_RATE"],
            "observed": round(map_rate, 4),
            "threshold": min_mapping_rate
        })

    status = "FAIL" if failed_gates else "PASS"
    return status, failed_gates, metrics


def main():
    parser = argparse.ArgumentParser(description="QC gate and no-call decision engine")
    parser.add_argument("--sample-id",         required=True)
    parser.add_argument("--fastp-json",         required=True)
    parser.add_argument("--bam-stats",          required=True)
    parser.add_argument("--ff-metrics",         required=True)
    parser.add_argument("--aneuploidy-calls",   required=True)
    parser.add_argument("--sca-calls",          required=True)
    parser.add_argument("--min-raw-reads",      type=int,   default=8000000)
    parser.add_argument("--min-unique-reads",   type=int,   default=4000000)
    parser.add_argument("--max-dup-rate",       type=float, default=0.35)
    parser.add_argument("--min-ff",             type=float, default=0.04)
    parser.add_argument("--min-mapping-rate",   type=float, default=0.70)
    parser.add_argument("--out",                required=True)
    args = parser.parse_args()

    fastp_metrics = parse_fastp_json(args.fastp_json)
    bam_metrics   = parse_bam_stats(args.bam_stats)
    ff_data       = parse_ff_metrics(args.ff_metrics)

    with open(args.aneuploidy_calls) as f:
        aneuploidy_data = json.load(f)
    with open(args.sca_calls) as f:
        sca_data = json.load(f)

    qc_status, failed_gates, qc_metrics = evaluate_qc_gates(
        fastp_metrics, bam_metrics, ff_data,
        args.min_raw_reads, args.min_unique_reads,
        args.max_dup_rate, args.min_ff, args.min_mapping_rate
    )

    # Build no-call reasons
    no_call_reasons = [g["code"] for g in failed_gates]

    # Determine reportability per target
    reportable = qc_status == "PASS"

    result = {
        "sample_id":       args.sample_id,
        "timestamp":       datetime.now(timezone.utc).isoformat(),
        "qc_status":       qc_status,
        "reportable":      reportable,
        "no_call_reasons": no_call_reasons,
        "failed_gates":    failed_gates,
        "qc_metrics":      qc_metrics,
        "thresholds": {
            "min_raw_reads":    args.min_raw_reads,
            "min_unique_reads": args.min_unique_reads,
            "max_dup_rate":     args.max_dup_rate,
            "min_ff":           args.min_ff,
            "min_mapping_rate": args.min_mapping_rate
        }
    }

    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    print(f"QC gate evaluation for {args.sample_id}: {qc_status}", file=sys.stderr)
    if failed_gates:
        for g in failed_gates:
            print(f"  FAILED [{g['code']}]: {g['message']} "
                  f"(observed={g['observed']}, threshold={g['threshold']})",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
