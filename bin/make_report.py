#!/usr/bin/env python3
"""
make_report.py
==============
EN: Generate clinical NIPT report:
    - Machine-readable JSON (validated against report schema)
    - Human-readable HTML summary
    - Audit trail JSON with MD5 checksums and environment metadata
VI: Tạo báo cáo lâm sàng NIPT:
    - JSON có thể đọc bằng máy (kiểm tra dựa trên schema báo cáo)
    - Tóm tắt HTML dành cho người dùng
    - JSON ghi lại dấu kiểm toán kèm MD5 checksum và thông tin môi trường
"""

import argparse
import json
import sys
import hashlib
import os
from datetime import datetime, timezone
from pathlib import Path


# EN: Clinical disclaimer — must appear on every report
# VI: Tuyên bố từ chối trách nhiệm lâm sàng — phải xuất hiện trên mọi báo cáo
CLINICAL_DISCLAIMER = (
    "NIPT is a screening test and is not diagnostic. A positive (high-risk) result "
    "does not confirm a chromosomal abnormality, and a negative (low-risk) result "
    "does not exclude one. All high-risk results should be confirmed by diagnostic "
    "testing (e.g., amniocentesis or chorionic villus sampling) following genetic "
    "counseling. / "
    "NIPT là xét nghiệm sàng lọc và KHÔNG phải xét nghiệm chẩn đoán. Kết quả cao "
    "nguy cơ không xác nhận bất thường NST, và kết quả thấp nguy cơ không loại trừ "
    "bất thường. Tất cả kết quả cao nguy cơ cần được xác nhận bằng xét nghiệm chẩn "
    "đoán (ví dụ: chọc ối hoặc sinh thiết gai nhau) sau tư vấn di truyền."
)

# EN: Display labels for classification codes
# VI: Nhãn hiển thị cho các mã phân loại kết quả
RISK_DISPLAY = {
    "HIGH_RISK":             "HIGH RISK / CAO NGUY CƠ",
    "LOW_RISK":              "LOW RISK / THẤP NGUY CƠ",
    "GREY_ZONE":             "BORDERLINE / VÙNG XÁM",
    "NO_CALL":               "NO CALL / KHÔNG GỌI ĐƯỢC",
    "SUSPECTED_RAA":         "SUSPECTED RAA / NGHI NGỜ RAA",
    "INSUFFICIENT_EVIDENCE": "INSUFFICIENT EVIDENCE / KHÔNG ĐỦ BẰNG CHỨNG"
}


def load_json(path: str) -> dict:
    """
    EN: Load a JSON file and return its content as a dict.
    VI: Tải file JSON và trả về nội dung dạng từ điển.
    """
    with open(path) as f:
        return json.load(f)


def file_md5(path: str) -> str:
    """
    EN: Compute MD5 hash of a file for audit trail integrity verification.
    VI: Tính MD5 hash của file để xác minh tính toàn vẹn trong dấu kiểm toán.
    """
    h = hashlib.md5()
    try:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()
    except (FileNotFoundError, IOError):
        return "N/A"


def build_report(args, qc: dict, aneuploidy: dict, sca: dict) -> dict:
    """
    EN: Build the structured clinical report dict from all module outputs.
    VI: Xây dựng từ điển báo cáo lâm sàng có cấu trúc từ tất cả đầu ra module.
    """
    now_iso = datetime.now(timezone.utc).isoformat()

    # EN: Aggregate per-target results
    # VI: Tổng hợp kết quả từng mục tiêu
    targets = {}

    # ── T21 / T18 / T13 ──────────────────────────────────────
    for chr_key, display_name in [("chr21", "T21"), ("chr18", "T18"), ("chr13", "T13")]:
        call = aneuploidy.get("targets", {}).get(chr_key, {})
        if not qc.get("reportable", False):
            # EN: Override to NO_CALL when QC gate failed
            # VI: Ghi đè thành NO_CALL khi cổng QC không đạt
            classification = "NO_CALL"
            reasons = qc.get("no_call_reasons", [])
        else:
            classification = call.get("classification", "NO_CALL")
            reasons = [call.get("no_call_reason")] if call.get("no_call_reason") else []

        targets[display_name] = {
            "chromosome":      chr_key,
            "classification":  classification,
            "clinical_risk":   RISK_DISPLAY.get(classification, classification),
            "z_score":         call.get("z_robust"),
            "no_call_reasons": [r for r in reasons if r]
        }

    # ── SCA chrX / chrY ───────────────────────────────────────
    for chr_key, display_name in [("chrX", "SCA_X"), ("chrY", "SCA_Y")]:
        call = sca.get("calls", {}).get(chr_key, {})
        if not qc.get("reportable", False):
            classification = "NO_CALL"
            reasons = qc.get("no_call_reasons", [])
        else:
            classification = call.get("classification", "NO_CALL")
            reasons = [call.get("no_call_reason")] if call.get("no_call_reason") else []

        targets[display_name] = {
            "chromosome":      chr_key,
            "classification":  classification,
            "clinical_risk":   RISK_DISPLAY.get(classification, classification),
            "z_score":         call.get("z"),
            "no_call_reasons": [r for r in reasons if r]
        }

    report = {
        "schema_version":    "1.0",
        "pipeline_version":  args.pipeline_version,
        "genome_build":      args.genome_build,
        "clinical_tier":     args.clinical_tier,
        "config_version":    args.config_version,
        "sample_id":         args.sample_id,
        "run_id":            args.run_id,
        "report_generated":  now_iso,
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
        "results":           targets,
        "inferred_sex":      sca.get("inferred_sex"),
        "reported_sex":      sca.get("reported_sex"),
        "disclaimer":        CLINICAL_DISCLAIMER
    }

    return report


def build_html(report: dict) -> str:
    """
    EN: Generate a simple HTML clinical report with color-coded risk levels.
    VI: Tạo báo cáo lâm sàng HTML đơn giản với các mức nguy cơ được tô màu.
    """
    def risk_class(risk: str) -> str:
        """EN: Map risk label to CSS class. | VI: Ánh xạ nhãn nguy cơ thành lớp CSS."""
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
    ff_val  = qc.get("fetal_fraction")
    ff_str  = f"{ff_val:.1%}" if ff_val is not None else "N/A"
    map_val = qc.get("mapping_rate")
    map_str = f"{map_val:.1%}" if map_val is not None else "N/A"

    # EN: Format integer counts safely — 'N/A' when missing/non-numeric
    # VI: Định dạng số đọc an toàn — hiển thị 'N/A' khi thiếu dữ liệu
    raw_reads_val    = qc.get("raw_reads")
    mapped_reads_val = qc.get("mapped_reads")
    raw_reads_str    = f"{raw_reads_val:,}" if isinstance(raw_reads_val, (int, float)) else "N/A"
    mapped_reads_str = f"{mapped_reads_val:,}" if isinstance(mapped_reads_val, (int, float)) else "N/A"

    html = f"""<!DOCTYPE html>
<html lang="vi-en">
<head>
<meta charset="UTF-8">
<!-- EN: Clinical NIPT report | VI: Báo cáo lâm sàng NIPT -->
<title>NIPT Clinical Report / Báo cáo lâm sàng NIPT - {report['sample_id']}</title>
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
<!-- EN: Report title | VI: Tiêu đề báo cáo -->
<h1>NIPT Clinical Report / Báo cáo Lâm Sàng NIPT</h1>

<!-- EN: Sample information section | VI: Phần thông tin mẫu -->
<h2>Sample Information / Thông tin Mẫu</h2>
<table>
  <tr><th>Field / Trường</th><th>Value / Giá trị</th></tr>
  <tr><td>Sample ID / Mã mẫu</td><td>{report['sample_id']}</td></tr>
  <tr><td>Run ID / Mã lần chạy</td><td>{report['run_id']}</td></tr>
  <tr><td>Report Generated / Ngày tạo báo cáo</td><td>{report['report_generated']}</td></tr>
  <tr><td>Pipeline Version / Phiên bản pipeline</td><td>{report['pipeline_version']}</td></tr>
  <tr><td>Genome Build / Phiên bản hệ gen</td><td>{report['genome_build']}</td></tr>
  <tr><td>Clinical Tier / Cấp độ lâm sàng</td><td>{report['clinical_tier']}</td></tr>
  <tr><td>Inferred Fetal Sex / Giới tính thai nhi suy luận</td><td>{report.get('inferred_sex', 'N/A')}</td></tr>
</table>

<!-- EN: QC summary section | VI: Phần tóm tắt QC -->
<h2>QC Summary / Tóm tắt Chất lượng</h2>
<table>
  <tr><th>Metric / Chỉ số</th><th>Value / Giá trị</th></tr>
  <tr><td>Overall QC / QC tổng thể</td><td class="{'qc-pass' if qc['status']=='PASS' else 'qc-fail'}">{qc['status']}</td></tr>
  <tr><td>Fetal Fraction / Tỷ lệ DNA thai nhi</td><td>{ff_str}</td></tr>
  <tr><td>Mapping Rate / Tỷ lệ ánh xạ</td><td>{map_str}</td></tr>
  <tr><td>Raw Reads / Số đọc thô</td><td>{raw_reads_str}</td></tr>
  <tr><td>Mapped Reads / Số đọc ánh xạ</td><td>{mapped_reads_str}</td></tr>
</table>

<!-- EN: Clinical results section | VI: Phần kết quả lâm sàng -->
<h2>Clinical Results / Kết quả Lâm sàng</h2>
<table>
  <tr>
    <th>Target / Mục tiêu</th>
    <th>Chromosome / NST</th>
    <th>Clinical Classification / Phân loại lâm sàng</th>
    <th>Z-Score</th>
    <th>Notes / Ghi chú</th>
  </tr>
  {rows}
</table>

<!-- EN: Clinical disclaimer | VI: Tuyên bố từ chối trách nhiệm lâm sàng -->
<div class="disclaimer">
  <strong>CLINICAL DISCLAIMER / TUYÊN BỐ LÂM SÀNG:</strong> {report['disclaimer']}
</div>

</body>
</html>"""
    return html


def main():
    # EN: Parse command-line arguments
    # VI: Phân tích đối số dòng lệnh
    parser = argparse.ArgumentParser(
        description="EN: Generate NIPT clinical report | VI: Tạo báo cáo lâm sàng NIPT"
    )
    parser.add_argument("--sample-id",        required=True,
                        help="EN: Sample identifier         | VI: Mã mẫu")
    parser.add_argument("--run-id",           required=True,
                        help="EN: Run identifier            | VI: Mã lần chạy")
    parser.add_argument("--qc-decision",      required=True,
                        help="EN: QC decision JSON          | VI: File JSON quyết định QC")
    parser.add_argument("--aneuploidy-calls", required=True,
                        help="EN: Aneuploidy calls JSON     | VI: File JSON kết quả lệch bội")
    parser.add_argument("--sca-calls",        required=True,
                        help="EN: SCA calls JSON            | VI: File JSON kết quả SCA")
    parser.add_argument("--pipeline-version", required=True,
                        help="EN: Pipeline version string   | VI: Chuỗi phiên bản pipeline")
    parser.add_argument("--genome-build",     required=True,
                        help="EN: Reference genome build    | VI: Phiên bản hệ gen tham chiếu")
    parser.add_argument("--clinical-tier",    required=True,
                        help="EN: Clinical tier (tier1/tier2) | VI: Cấp độ lâm sàng")
    parser.add_argument("--config-version",   required=True,
                        help="EN: Config version string     | VI: Phiên bản cấu hình")
    parser.add_argument("--report-schema",    required=True,
                        help="EN: Report JSON schema file   | VI: File JSON schema báo cáo")
    parser.add_argument("--out-json",         required=True,
                        help="EN: Output report JSON        | VI: File JSON báo cáo đầu ra")
    parser.add_argument("--out-html",         required=True,
                        help="EN: Output report HTML        | VI: File HTML báo cáo đầu ra")
    parser.add_argument("--out-audit",        required=True,
                        help="EN: Output audit trail JSON   | VI: File JSON dấu kiểm toán đầu ra")
    args = parser.parse_args()

    print(f"[INFO] EN: Generating clinical report / VI: Đang tạo báo cáo lâm sàng "
          f"cho mẫu: {args.sample_id} ...", file=sys.stderr)

    # ── Load all module outputs / Tải tất cả đầu ra module ───
    qc         = load_json(args.qc_decision)
    aneuploidy = load_json(args.aneuploidy_calls)
    sca        = load_json(args.sca_calls)

    # ── Build clinical report / Xây dựng báo cáo lâm sàng ───
    report = build_report(args, qc, aneuploidy, sca)

    # EN: Optionally validate against JSON schema
    # VI: Tùy chọn kiểm tra báo cáo dựa trên JSON schema
    try:
        import jsonschema
        with open(args.report_schema) as f:
            schema = json.load(f)
        jsonschema.validate(instance=report, schema=schema)
        print(f"[INFO] EN: Report schema validation passed / "
              f"VI: Kiểm tra schema báo cáo đạt", file=sys.stderr)
    except ImportError:
        pass  # EN: jsonschema not installed — skip | VI: chưa cài jsonschema — bỏ qua
    except jsonschema.ValidationError as e:
        print(f"WARNING / CẢNH BÁO: Report schema validation failed / "
              f"Kiểm tra schema báo cáo thất bại: {e.message}", file=sys.stderr)

    # ── Write JSON report / Ghi báo cáo JSON ─────────────────
    with open(args.out_json, "w") as f:
        json.dump(report, f, indent=2)

    # ── Write HTML report (intentional clinical output file)
    # ── Ghi báo cáo HTML (file đầu ra lâm sàng có chủ đích)
    html_content = build_html(report)
    with open(args.out_html, "w") as f:
        f.write(html_content)

    # ── Write audit trail / Ghi dấu kiểm toán ────────────────
    # EN: Audit trail records all inputs, versions, and environment metadata
    # VI: Dấu kiểm toán ghi lại tất cả đầu vào, phiên bản và siêu dữ liệu môi trường
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

    print(f"[DONE / HOÀN THÀNH] EN: Report generated / VI: Đã tạo báo cáo "
          f"cho mẫu {args.sample_id}: QC={qc.get('qc_status', 'UNKNOWN')}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
