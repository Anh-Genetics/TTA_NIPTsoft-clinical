#!/usr/bin/env python3
"""
qc_decision.py
==============
EN: Unified QC gate and no-call decision engine.
VI: Bộ máy kiểm tra QC thống nhất và quyết định "không gọi được kết quả".

EN: Evaluates all QC metrics and determines:
    1. Overall QC status (PASS / FAIL)
    2. No-call reasons (if any)
    3. Whether each target result should be reported or suppressed
VI: Đánh giá tất cả chỉ số QC và xác định:
    1. Trạng thái QC tổng thể (ĐẠT / KHÔNG ĐẠT)
    2. Lý do không gọi được kết quả (nếu có)
    3. Liệu kết quả từng mục tiêu có được báo cáo hay không

Decision gates / Các cổng kiểm tra:
- G01: Minimum raw reads      / Số đọc thô tối thiểu
- G02: Minimum unique reads   / Số đọc ánh xạ duy nhất tối thiểu
- G03: Maximum duplication    / Tỷ lệ trùng lặp tối đa
- G04: Minimum fetal fraction / Tỷ lệ DNA thai nhi tối thiểu
- G05: Minimum mapping rate   / Tỷ lệ ánh xạ tối thiểu
- G06: Missing fetal fraction / Tỷ lệ DNA thai nhi không tính được
"""

import argparse
import json
import sys
from datetime import datetime, timezone


# EN: Gate codes with English descriptions (Vietnamese in parentheses)
# VI: Mã cổng kiểm tra kèm mô tả tiếng Anh (tiếng Việt trong ngoặc)
QC_GATE_CODES = {
    "G01_LOW_RAW_READS":     "Raw read count below minimum threshold "
                             "(Số đọc thô dưới ngưỡng tối thiểu)",
    "G02_LOW_UNIQUE_READS":  "Unique mapped reads below minimum threshold "
                             "(Số đọc ánh xạ duy nhất dưới ngưỡng tối thiểu)",
    "G03_HIGH_DUP_RATE":     "Duplication rate exceeds maximum threshold "
                             "(Tỷ lệ trùng lặp vượt quá ngưỡng tối đa)",
    "G04_LOW_FF":            "Fetal fraction below minimum threshold "
                             "(Tỷ lệ DNA thai nhi dưới ngưỡng tối thiểu)",
    "G05_LOW_MAPPING_RATE":  "Mapping rate below minimum threshold "
                             "(Tỷ lệ ánh xạ dưới ngưỡng tối thiểu)",
    "G06_MISSING_FF":        "Fetal fraction could not be estimated "
                             "(Không thể ước tính tỷ lệ DNA thai nhi)",
    "G07_MISSING_BAM_STATS": "BAM statistics could not be parsed "
                             "(Không thể phân tích số liệu thống kê BAM)"
}


def parse_fastp_json(path: str) -> dict:
    """
    EN: Extract summary metrics from fastp JSON output.
    VI: Trích xuất các chỉ số tóm tắt từ file JSON của fastp.
    """
    try:
        with open(path) as f:
            data = json.load(f)
        summary = data.get("summary", {})
        before  = summary.get("before_filtering", {})
        after   = summary.get("after_filtering",  {})
        return {
            "raw_reads":        before.get("total_reads", 0),
            "clean_reads":      after.get("total_reads",  0),
            "q20_rate":         after.get("q20_rate",     0.0),
            "q30_rate":         after.get("q30_rate",     0.0),
            "gc_content":       after.get("gc_content",   0.0),
            "duplication_rate": data.get("duplication", {}).get("rate", 0.0)
        }
    except (FileNotFoundError, json.JSONDecodeError, KeyError):
        # EN: Return zeros so downstream gates can still evaluate
        # VI: Trả về giá trị zero để các cổng tiếp theo vẫn có thể đánh giá
        return {
            "raw_reads": 0, "clean_reads": 0,
            "q20_rate":  0.0, "q30_rate":  0.0,
            "gc_content": 0.0, "duplication_rate": 0.0
        }


def parse_bam_stats(path: str) -> dict:
    """
    EN: Parse samtools flagstat output into a metrics dict.
    VI: Phân tích file đầu ra của samtools flagstat thành từ điển chỉ số.
    """
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
        # EN: Derive rates from raw counts
        # VI: Tính tỷ lệ từ số đếm thô
        if metrics["total_reads"] > 0:
            metrics["mapping_rate"] = metrics["mapped_reads"] / metrics["total_reads"]
        if metrics["mapped_reads"] > 0:
            metrics["dup_rate"] = metrics["duplicate_reads"] / metrics["mapped_reads"]
    except (FileNotFoundError, ValueError, IndexError):
        pass
    return metrics


def parse_ff_metrics(path: str) -> dict:
    """
    EN: Parse fetal fraction metrics JSON.
    VI: Phân tích file JSON chứa các chỉ số tỷ lệ DNA thai nhi.
    """
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
    EN: Evaluate all QC gates. Returns (status, failed_gates, all_metrics).
    VI: Đánh giá tất cả cổng QC. Trả về (trạng thái, cổng thất bại, tất cả chỉ số).
    """
    failed_gates = []

    # EN: Extract individual metrics from parsed dicts
    # VI: Trích xuất từng chỉ số từ các từ điển đã phân tích
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

    # ── Gate G01: Raw reads / Số đọc thô ─────────────────────
    if raw_reads < min_raw_reads:
        failed_gates.append({
            "code":      "G01_LOW_RAW_READS",
            "message":   QC_GATE_CODES["G01_LOW_RAW_READS"],
            "observed":  raw_reads,
            "threshold": min_raw_reads
        })

    # ── Gate G02: Unique mapped reads / Số đọc ánh xạ duy nhất ─
    if mapped_reads < min_unique_reads:
        failed_gates.append({
            "code":      "G02_LOW_UNIQUE_READS",
            "message":   QC_GATE_CODES["G02_LOW_UNIQUE_READS"],
            "observed":  mapped_reads,
            "threshold": min_unique_reads
        })

    # ── Gate G03: Duplication rate / Tỷ lệ trùng lặp ─────────
    if dup_rate > max_dup_rate:
        failed_gates.append({
            "code":      "G03_HIGH_DUP_RATE",
            "message":   QC_GATE_CODES["G03_HIGH_DUP_RATE"],
            "observed":  round(dup_rate, 4),
            "threshold": max_dup_rate
        })

    # ── Gate G04 / G06: Fetal fraction / Tỷ lệ DNA thai nhi ──
    if ff_val is None:
        failed_gates.append({
            "code":      "G06_MISSING_FF",
            "message":   QC_GATE_CODES["G06_MISSING_FF"],
            "observed":  None,
            "threshold": min_ff
        })
    elif ff_val < min_ff:
        failed_gates.append({
            "code":      "G04_LOW_FF",
            "message":   QC_GATE_CODES["G04_LOW_FF"],
            "observed":  round(ff_val, 4),
            "threshold": min_ff
        })

    # ── Gate G05: Mapping rate / Tỷ lệ ánh xạ ────────────────
    # EN: Only evaluate when we have actual BAM data (mapped_reads > 0)
    # VI: Chỉ đánh giá khi có dữ liệu BAM thực tế (mapped_reads > 0)
    if map_rate < min_mapping_rate and mapped_reads > 0:
        failed_gates.append({
            "code":      "G05_LOW_MAPPING_RATE",
            "message":   QC_GATE_CODES["G05_LOW_MAPPING_RATE"],
            "observed":  round(map_rate, 4),
            "threshold": min_mapping_rate
        })

    status = "FAIL" if failed_gates else "PASS"
    return status, failed_gates, metrics


def main():
    # EN: Parse command-line arguments
    # VI: Phân tích đối số dòng lệnh
    parser = argparse.ArgumentParser(
        description="EN: QC gate and no-call decision engine | "
                    "VI: Bộ máy kiểm tra QC và quyết định không gọi kết quả"
    )
    parser.add_argument("--sample-id",         required=True,
                        help="EN: Sample identifier | VI: Mã mẫu")
    parser.add_argument("--fastp-json",         required=True,
                        help="EN: fastp JSON output | VI: File JSON đầu ra của fastp")
    parser.add_argument("--bam-stats",          required=True,
                        help="EN: samtools flagstat output | VI: File đầu ra flagstat")
    parser.add_argument("--ff-metrics",         required=True,
                        help="EN: Fetal fraction JSON | VI: File JSON tỷ lệ DNA thai nhi")
    parser.add_argument("--aneuploidy-calls",   required=True,
                        help="EN: Aneuploidy calls JSON | VI: File JSON kết quả gọi lệch bội")
    parser.add_argument("--sca-calls",          required=True,
                        help="EN: SCA calls JSON | VI: File JSON kết quả gọi lệch bội NST giới tính")
    parser.add_argument("--min-raw-reads",      type=int,   default=8000000,
                        help="EN: Min raw reads (default 8M) | VI: Số đọc thô tối thiểu (mặc định 8 triệu)")
    parser.add_argument("--min-unique-reads",   type=int,   default=4000000,
                        help="EN: Min unique mapped reads (default 4M) | VI: Số đọc ánh xạ duy nhất tối thiểu")
    parser.add_argument("--max-dup-rate",       type=float, default=0.35,
                        help="EN: Max duplication rate (default 0.35) | VI: Tỷ lệ trùng lặp tối đa (mặc định 0.35)")
    parser.add_argument("--min-ff",             type=float, default=0.04,
                        help="EN: Min fetal fraction (default 0.04) | VI: Tỷ lệ DNA thai nhi tối thiểu (mặc định 4%)")
    parser.add_argument("--min-mapping-rate",   type=float, default=0.70,
                        help="EN: Min mapping rate (default 0.70) | VI: Tỷ lệ ánh xạ tối thiểu (mặc định 70%)")
    parser.add_argument("--out",                required=True,
                        help="EN: Output JSON file | VI: File JSON kết quả đầu ra")
    args = parser.parse_args()

    # ── Load all inputs / Tải tất cả dữ liệu đầu vào ─────────
    print(f"[INFO] EN: Evaluating QC gates / VI: Đang đánh giá cổng QC "
          f"cho mẫu: {args.sample_id} ...", file=sys.stderr)

    fastp_metrics = parse_fastp_json(args.fastp_json)
    bam_metrics   = parse_bam_stats(args.bam_stats)
    ff_data       = parse_ff_metrics(args.ff_metrics)

    with open(args.aneuploidy_calls) as f:
        aneuploidy_data = json.load(f)
    with open(args.sca_calls) as f:
        sca_data = json.load(f)

    # ── Run QC evaluation / Chạy đánh giá QC ─────────────────
    qc_status, failed_gates, qc_metrics = evaluate_qc_gates(
        fastp_metrics, bam_metrics, ff_data,
        args.min_raw_reads, args.min_unique_reads,
        args.max_dup_rate, args.min_ff, args.min_mapping_rate
    )

    # EN: Build no-call reason codes from failed gates
    # VI: Xây dựng danh sách mã lý do không gọi được từ các cổng thất bại
    no_call_reasons = [g["code"] for g in failed_gates]

    # EN: Sample is reportable only when ALL gates pass
    # VI: Mẫu chỉ được báo cáo khi TẤT CẢ cổng đều đạt
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

    # ── Write output / Ghi kết quả ────────────────────────────
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    # EN: Print QC summary to stderr
    # VI: In tóm tắt QC ra stderr
    print(f"[QC] EN: Gate evaluation result / VI: Kết quả đánh giá cổng QC "
          f"cho mẫu {args.sample_id}: {qc_status}", file=sys.stderr)
    if failed_gates:
        for g in failed_gates:
            print(f"  [FAIL] [{g['code']}]: {g['message']} "
                  f"(observed={g['observed']}, threshold={g['threshold']})",
                  file=sys.stderr)

    print(f"[DONE / HOÀN THÀNH] QC decision / Quyết định QC "
          f"cho mẫu: {args.sample_id}", file=sys.stderr)


if __name__ == "__main__":
    main()
