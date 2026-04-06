#!/usr/bin/env python3
"""
make_report.py
Generate clinical NIPT report:
- Machine-readable JSON (validated against report schema)
- Human-readable HTML summary
- Audit trail JSON
"""

import argparse
import json
import sys
import hashlib
import os
from datetime import datetime, timezone
from pathlib import Path


CLINICAL_DISCLAIMER = (
    "NIPT is a screening test and is not diagnostic. A positive (high-risk) result "
    "does not confirm a chromosomal abnormality, and a negative (low-risk) result "
    "does not exclude one. All high-risk results should be confirmed by diagnostic "
    "testing (e.g., amniocentesis or chorionic villus sampling) following genetic "
    "counseling."
)

RISK_DISPLAY = {
    "HIGH_RISK":       "HIGH RISK",
    "LOW_RISK":        "LOW RISK",
    "GREY_ZONE":       "BORDERLINE",
    "NO_CALL":         "NO CALL",
    "SUSPECTED_RAA":   "SUSPECTED RAA",
    "INSUFFICIENT_EVIDENCE": "INSUFFICIENT EVIDENCE"
}


def load_json(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def file_md5(path: str) -> str:
    """Compute MD5 of a file for audit trail."""
    h = hashlib.md5()
    try:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()
    except (FileNotFoundError, IOError):
        return "N/A"


def build_report(args, qc: dict, aneuploidy: dict, sca: dict) -> dict:
    """Build the structured clinical report."""
    now_iso = datetime.now(timezone.utc).isoformat()

    # Resolve per-target results
    targets = {}

    # T21 / T18 / T13
    for chr_key, display_name in [("chr21", "T21"), ("chr18", "T18"), ("chr13", "T13")]:
        call = aneuploidy.get("targets", {}).get(chr_key, {})
        if not qc.get("reportable", False):
            classification = "NO_CALL"
            reasons = qc.get("no_call_reasons", [])
        else:
            classification = call.get("classification", "NO_CALL")
            reasons = [call.get("no_call_reason")] if call.get("no_call_reason") else []

        targets[display_name] = {
            "chromosome":     chr_key,
            "classification": classification,
            "clinical_risk":  RISK_DISPLAY.get(classification, classification),
            "z_score":        call.get("z_robust"),
            "no_call_reasons": [r for r in reasons if r]
        }

    # SCA (chrX / chrY)
    for chr_key, display_name in [("chrX", "SCA_X"), ("chrY", "SCA_Y")]:
        call = sca.get("calls", {}).get(chr_key, {})
        if not qc.get("reportable", False):
            classification = "NO_CALL"
            reasons = qc.get("no_call_reasons", [])
        else:
            classification = call.get("classification", "NO_CALL")
            reasons = [call.get("no_call_reason")] if call.get("no_call_reason") else []

        targets[display_name] = {
            "chromosome":     chr_key,
            "classification": classification,
            "clinical_risk":  RISK_DISPLAY.get(classification, classification),
            "z_score":        call.get("z"),
            "no_call_reasons": [r for r in reasons if r]
        }

    report = {
        "schema_version":   "1.0",
        "pipeline_version": args.pipeline_version,
        "genome_build":     args.genome_build,
        "clinical_tier":    args.clinical_tier,
        "config_version":   args.config_version,
        "sample_id":        args.sample_id,
        "run_id":           args.run_id,
        "report_generated": now_iso,
        "qc_summary": {
            "status":          qc.get("qc_status", "UNKNOWN"),
            "reportable":      qc.get("reportable", False),
            "fetal_fraction":  qc.get("qc_metrics", {}).get("fetal_fraction"),
            "mapping_rate":    qc.get("qc_metrics", {}).get("mapping_rate"),
            "raw_reads":       qc.get("qc_metrics", {}).get("raw_reads"),
            "mapped_reads":    qc.get("qc_metrics", {}).get("mapped_reads"),
            "dup_rate":        qc.get("qc_metrics", {}).get("dup_rate"),
            "no_call_reasons": qc.get("no_call_reasons", []),
            "failed_gates":    qc.get("failed_gates", [])
        },
        "results":          targets,
        "inferred_sex":     sca.get("inferred_sex"),
        "reported_sex":     sca.get("reported_sex"),
        "disclaimer":       CLINICAL_DISCLAIMER
    }

    return report


def build_html(report: dict) -> str:
    """Generate simple HTML report."""
    def risk_class(risk: str) -> str:
        if "HIGH" in risk:    return "high-risk"
        if "LOW" in risk:     return "low-risk"
        if "NO CALL" in risk: return "no-call"
        return "borderline"

    rows = ""
    for target, info in report["results"].items():
        clinical_risk = info["clinical_risk"]
        z = f"{info['z_score']:.2f}" if info.get("z_score") is not None else "N/A"
        reasons = "; ".join(info.get("no_call_reasons", [])) or ""
        rows += f"""
        <tr>
            <td>{target}</td>
            <td>{info['chromosome']}</td>
            <td class="{risk_class(clinical_risk)}">{clinical_risk}</td>
            <td>{z}</td>
            <td>{reasons}</td>
        </tr>"""

    qc = report["qc_summary"]
    ff_val = qc.get("fetal_fraction")
    ff_str = f"{ff_val:.1%}" if ff_val is not None else "N/A"
    map_val = qc.get("mapping_rate")
    map_str = f"{map_val:.1%}" if map_val is not None else "N/A"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>NIPT Clinical Report - {report['sample_id']}</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; color: #222; }}
  h1   {{ color: #003366; }}
  h2   {{ color: #005599; border-bottom: 1px solid #ccc; padding-bottom: 5px; }}
  table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
  th   {{ background: #003366; color: white; padding: 8px 12px; text-align: left; }}
  td   {{ padding: 7px 12px; border-bottom: 1px solid #eee; }}
  tr:hover {{ background: #f5f5f5; }}
  .high-risk  {{ color: #c00; font-weight: bold; }}
  .low-risk   {{ color: #060; font-weight: bold; }}
  .no-call    {{ color: #888; font-style: italic; }}
  .borderline {{ color: #c80; font-weight: bold; }}
  .disclaimer {{ background: #fff8e1; border-left: 4px solid #f9a825;
                 padding: 12px; margin-top: 20px; font-size: 0.9em; }}
  .qc-pass {{ color: #060; font-weight: bold; }}
  .qc-fail {{ color: #c00; font-weight: bold; }}
</style>
</head>
<body>
<h1>NIPT Clinical Report</h1>

<h2>Sample Information</h2>
<table>
  <tr><th>Field</th><th>Value</th></tr>
  <tr><td>Sample ID</td><td>{report['sample_id']}</td></tr>
  <tr><td>Run ID</td><td>{report['run_id']}</td></tr>
  <tr><td>Report Generated</td><td>{report['report_generated']}</td></tr>
  <tr><td>Pipeline Version</td><td>{report['pipeline_version']}</td></tr>
  <tr><td>Genome Build</td><td>{report['genome_build']}</td></tr>
  <tr><td>Clinical Tier</td><td>{report['clinical_tier']}</td></tr>
  <tr><td>Inferred Fetal Sex</td><td>{report.get('inferred_sex', 'N/A')}</td></tr>
</table>

<h2>QC Summary</h2>
<table>
  <tr><th>Metric</th><th>Value</th></tr>
  <tr><td>Overall QC</td><td class="{'qc-pass' if qc['status']=='PASS' else 'qc-fail'}">{qc['status']}</td></tr>
  <tr><td>Fetal Fraction</td><td>{ff_str}</td></tr>
  <tr><td>Mapping Rate</td><td>{map_str}</td></tr>
  <tr><td>Raw Reads</td><td>{qc.get('raw_reads', 'N/A'):,}</td></tr>
  <tr><td>Mapped Reads</td><td>{qc.get('mapped_reads', 'N/A'):,}</td></tr>
</table>

<h2>Clinical Results</h2>
<table>
  <tr>
    <th>Target</th>
    <th>Chromosome</th>
    <th>Clinical Classification</th>
    <th>Z-Score</th>
    <th>Notes</th>
  </tr>
  {rows}
</table>

<div class="disclaimer">
  <strong>CLINICAL DISCLAIMER:</strong> {report['disclaimer']}
</div>

</body>
</html>"""
    return html


def main():
    parser = argparse.ArgumentParser(description="Generate NIPT clinical report")
    parser.add_argument("--sample-id",        required=True)
    parser.add_argument("--run-id",           required=True)
    parser.add_argument("--qc-decision",      required=True)
    parser.add_argument("--aneuploidy-calls", required=True)
    parser.add_argument("--sca-calls",        required=True)
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--genome-build",     required=True)
    parser.add_argument("--clinical-tier",    required=True)
    parser.add_argument("--config-version",   required=True)
    parser.add_argument("--report-schema",    required=True)
    parser.add_argument("--out-json",         required=True)
    parser.add_argument("--out-html",         required=True)
    parser.add_argument("--out-audit",        required=True)
    args = parser.parse_args()

    qc         = load_json(args.qc_decision)
    aneuploidy = load_json(args.aneuploidy_calls)
    sca        = load_json(args.sca_calls)

    report = build_report(args, qc, aneuploidy, sca)

    # Optional: validate against schema
    try:
        import jsonschema
        with open(args.report_schema) as f:
            schema = json.load(f)
        jsonschema.validate(instance=report, schema=schema)
    except ImportError:
        pass  # jsonschema not available
    except jsonschema.ValidationError as e:
        print(f"WARNING: Report schema validation failed: {e.message}", file=sys.stderr)

    with open(args.out_json, "w") as f:
        json.dump(report, f, indent=2)

    # Write HTML report — intentional clinical output file
    html_content = build_html(report)
    with open(args.out_html, "w") as f:
        f.write(html_content)

    # Audit trail
    audit = {
        "sample_id":        args.sample_id,
        "run_id":           args.run_id,
        "timestamp":        datetime.now(timezone.utc).isoformat(),
        "pipeline_version": args.pipeline_version,
        "genome_build":     args.genome_build,
        "clinical_tier":    args.clinical_tier,
        "config_version":   args.config_version,
        "input_files": {
            "qc_decision":      {"path": args.qc_decision,      "md5": file_md5(args.qc_decision)},
            "aneuploidy_calls": {"path": args.aneuploidy_calls, "md5": file_md5(args.aneuploidy_calls)},
            "sca_calls":        {"path": args.sca_calls,        "md5": file_md5(args.sca_calls)}
        },
        "output_files": {
            "report_json": args.out_json,
            "report_html": args.out_html
        },
        "environment": {
            "python_version": sys.version,
            "hostname":       os.uname().nodename
        }
    }

    with open(args.out_audit, "w") as f:
        json.dump(audit, f, indent=2)

    print(f"Report generated for {args.sample_id}: QC={qc.get('qc_status', 'UNKNOWN')}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
