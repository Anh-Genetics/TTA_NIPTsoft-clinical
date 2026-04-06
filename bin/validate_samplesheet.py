#!/usr/bin/env python3
"""
validate_samplesheet.py
=======================
EN: Validate a NIPT samplesheet CSV against a JSON schema.
    Outputs a cleaned, validated CSV and a validation report JSON.
VI: Kiểm tra file samplesheet NIPT (CSV) dựa trên JSON schema.
    Xuất file CSV đã được làm sạch và kiểm tra, cùng báo cáo kiểm tra JSON.
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from datetime import datetime, timezone


def load_schema(schema_path: str) -> dict:
    """
    EN: Load JSON schema from file.
    VI: Tải JSON schema từ file.
    """
    with open(schema_path) as f:
        return json.load(f)


def validate_row(row: dict, schema: dict, row_num: int, errors: list) -> dict:
    """
    EN: Validate a single samplesheet row against schema properties.
        Appends error dicts to the errors list in place.
    VI: Kiểm tra một hàng của samplesheet dựa trên các thuộc tính schema.
        Thêm các lỗi vào danh sách errors tại chỗ.
    """
    required_fields = schema.get("required", [])
    properties      = schema.get("properties", {})

    # EN: Check all required fields are present and non-empty
    # VI: Kiểm tra tất cả trường bắt buộc có mặt và không trống
    for field in required_fields:
        if field not in row or not row[field].strip():
            errors.append({
                "row":     row_num,
                "field":   field,
                "message": f"Required field '{field}' is missing or empty "
                           f"(Trường bắt buộc '{field}' bị thiếu hoặc trống)"
            })

    # EN: Validate type, enum constraints, and FASTQ file existence
    # VI: Kiểm tra kiểu dữ liệu, ràng buộc enum và sự tồn tại của file FASTQ
    for field, spec in properties.items():
        val = row.get(field, "").strip()
        if not val:
            continue

        field_type    = spec.get("type")
        field_enum    = spec.get("enum")
        field_pattern = spec.get("pattern")

        # EN: Enum constraint check
        # VI: Kiểm tra ràng buộc enum (danh sách giá trị cho phép)
        if field_enum and val not in field_enum:
            errors.append({
                "row":     row_num,
                "field":   field,
                "message": f"Value '{val}' not in allowed values / "
                           f"Giá trị '{val}' không nằm trong danh sách cho phép: {field_enum}"
            })

        # EN: Numeric type check
        # VI: Kiểm tra kiểu số
        if field_type == "number":
            try:
                float(val)
            except ValueError:
                errors.append({
                    "row":     row_num,
                    "field":   field,
                    "message": f"Expected numeric value / Cần giá trị số, got: '{val}'"
                })

        # EN: FASTQ file existence check (warns — file may be on remote)
        # VI: Kiểm tra sự tồn tại file FASTQ (cảnh báo — có thể ở hệ thống từ xa)
        if field == "fastq_path":
            if not Path(val).exists():
                errors.append({
                    "row":     row_num,
                    "field":   field,
                    "message": f"FASTQ file not found / Không tìm thấy file FASTQ: {val}"
                })

    # EN: Return cleaned row with whitespace stripped
    # VI: Trả về hàng đã được làm sạch (loại bỏ khoảng trắng thừa)
    return {k: v.strip() if v else v for k, v in row.items()}


def main():
    # EN: Parse command-line arguments
    # VI: Phân tích đối số dòng lệnh
    parser = argparse.ArgumentParser(
        description="EN: Validate NIPT samplesheet | VI: Kiểm tra samplesheet NIPT"
    )
    parser.add_argument("--samplesheet",  required=True,
                        help="EN: Input CSV samplesheet | VI: File CSV samplesheet đầu vào")
    parser.add_argument("--schema",       required=True,
                        help="EN: JSON schema file      | VI: File JSON schema")
    parser.add_argument("--out-csv",      required=True,
                        help="EN: Output validated CSV  | VI: File CSV đã kiểm tra đầu ra")
    parser.add_argument("--out-report",   required=True,
                        help="EN: Output validation report JSON | VI: File JSON báo cáo kiểm tra")
    args = parser.parse_args()

    print(f"[INFO] EN: Validating samplesheet / VI: Đang kiểm tra samplesheet: "
          f"{args.samplesheet} ...", file=sys.stderr)

    schema = load_schema(args.schema)
    errors        = []
    warnings      = []
    validated_rows = []

    # EN: Read all rows from the CSV samplesheet
    # VI: Đọc tất cả hàng từ file CSV samplesheet
    with open(args.samplesheet) as f:
        reader = csv.DictReader(f)
        rows   = list(reader)

    if not rows:
        errors.append({"row": 0, "field": "file",
                       "message": "Samplesheet is empty / Samplesheet trống"})
    else:
        seen_ids = {}
        for i, row in enumerate(rows, start=2):  # EN: row 1 = header | VI: hàng 1 = tiêu đề
            cleaned = validate_row(row, schema, i, errors)

            # EN: Check for duplicate sample IDs
            # VI: Kiểm tra mã mẫu trùng lặp
            sid = cleaned.get("sample_id", "")
            if sid in seen_ids:
                errors.append({
                    "row":     i,
                    "field":   "sample_id",
                    "message": f"Duplicate sample_id / Mã mẫu trùng lặp: '{sid}' "
                               f"(first seen at row / lần đầu thấy ở hàng {seen_ids[sid]})"
                })
            else:
                seen_ids[sid] = i
            validated_rows.append(cleaned)

    # EN: Assemble validation report
    # VI: Tổng hợp báo cáo kiểm tra
    report = {
        "timestamp":      datetime.now(timezone.utc).isoformat(),
        "samplesheet":    args.samplesheet,
        "total_samples":  len(validated_rows),
        "errors":         errors,
        "warnings":       warnings,
        "status":         "FAIL" if errors else "PASS"
    }

    with open(args.out_report, "w") as f:
        json.dump(report, f, indent=2)

    # EN: Hard fail on errors — pipeline should not proceed with invalid samples
    # VI: Lỗi nghiêm trọng — pipeline không được tiếp tục với mẫu không hợp lệ
    if errors:
        print(f"ERROR / LỖI: Samplesheet validation FAILED / THẤT BẠI "
              f"with {len(errors)} error(s) / lỗi:", file=sys.stderr)
        for e in errors:
            print(f"  Row / Hàng {e['row']}: [{e['field']}] {e['message']}", file=sys.stderr)
        sys.exit(1)

    # EN: Write cleaned CSV when all rows pass
    # VI: Ghi CSV đã làm sạch khi tất cả hàng đều hợp lệ
    if validated_rows:
        with open(args.out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=validated_rows[0].keys())
            writer.writeheader()
            writer.writerows(validated_rows)
    else:
        Path(args.out_csv).touch()

    print(f"[DONE / HOÀN THÀNH] EN: Samplesheet validation PASSED / ĐẠT: "
          f"{len(validated_rows)} sample(s) / mẫu", file=sys.stderr)


if __name__ == "__main__":
    main()
