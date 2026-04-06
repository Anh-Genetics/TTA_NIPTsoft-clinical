#!/usr/bin/env python3
"""
validate_samplesheet.py
Validate a NIPT samplesheet CSV against a JSON schema.
Outputs a cleaned, validated CSV and a validation report JSON.
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from datetime import datetime, timezone


def load_schema(schema_path: str) -> dict:
    with open(schema_path) as f:
        return json.load(f)


def validate_row(row: dict, schema: dict, row_num: int, errors: list) -> dict:
    """Validate a single samplesheet row against schema properties."""
    required_fields = schema.get("required", [])
    properties = schema.get("properties", {})

    for field in required_fields:
        if field not in row or not row[field].strip():
            errors.append({
                "row": row_num,
                "field": field,
                "message": f"Required field '{field}' is missing or empty"
            })

    # Type/format validation
    for field, spec in properties.items():
        val = row.get(field, "").strip()
        if not val:
            continue

        field_type = spec.get("type")
        field_enum = spec.get("enum")
        field_pattern = spec.get("pattern")

        if field_enum and val not in field_enum:
            errors.append({
                "row": row_num,
                "field": field,
                "message": f"Value '{val}' not in allowed values: {field_enum}"
            })

        if field_type == "number":
            try:
                float(val)
            except ValueError:
                errors.append({
                    "row": row_num,
                    "field": field,
                    "message": f"Expected numeric value, got: '{val}'"
                })

        if field == "fastq_path":
            if not Path(val).exists():
                errors.append({
                    "row": row_num,
                    "field": field,
                    "message": f"FASTQ file not found: {val}"
                })

    cleaned = {}
    for k, v in row.items():
        cleaned[k] = v.strip() if v else v
    return cleaned


def main():
    parser = argparse.ArgumentParser(description="Validate NIPT samplesheet")
    parser.add_argument("--samplesheet",  required=True, help="Input CSV samplesheet")
    parser.add_argument("--schema",       required=True, help="JSON schema file")
    parser.add_argument("--out-csv",      required=True, help="Output validated CSV")
    parser.add_argument("--out-report",   required=True, help="Output validation report JSON")
    args = parser.parse_args()

    schema = load_schema(args.schema)
    errors = []
    warnings = []
    validated_rows = []

    with open(args.samplesheet) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        errors.append({"row": 0, "field": "file", "message": "Samplesheet is empty"})
    else:
        seen_ids = {}
        for i, row in enumerate(rows, start=2):  # row 1 = header
            cleaned = validate_row(row, schema, i, errors)
            # Check for duplicate sample_id
            sid = cleaned.get("sample_id", "")
            if sid in seen_ids:
                errors.append({
                    "row": i,
                    "field": "sample_id",
                    "message": f"Duplicate sample_id: '{sid}' (first seen at row {seen_ids[sid]})"
                })
            else:
                seen_ids[sid] = i
            validated_rows.append(cleaned)

    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "samplesheet": args.samplesheet,
        "total_samples": len(validated_rows),
        "errors": errors,
        "warnings": warnings,
        "status": "FAIL" if errors else "PASS"
    }

    with open(args.out_report, "w") as f:
        json.dump(report, f, indent=2)

    if errors:
        print(f"ERROR: Samplesheet validation failed with {len(errors)} error(s):", file=sys.stderr)
        for e in errors:
            print(f"  Row {e['row']}: [{e['field']}] {e['message']}", file=sys.stderr)
        sys.exit(1)

    # Write cleaned CSV
    if validated_rows:
        with open(args.out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=validated_rows[0].keys())
            writer.writeheader()
            writer.writerows(validated_rows)
    else:
        Path(args.out_csv).touch()

    print(f"Samplesheet validation PASSED: {len(validated_rows)} sample(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
