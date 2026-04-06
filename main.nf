#!/usr/bin/env nextflow
// ============================================================
// TTA_NIPTsoft-clinical  — Illumina single-end NIPT pipeline
// EN: Nextflow DSL2 main orchestrator — coordinates all pipeline stages
// VI: Bộ điều phối chính Nextflow DSL2 — phối hợp tất cả giai đoạn pipeline
// ============================================================
nextflow.enable.dsl = 2

// ------------------------------------------------------------
// EN: Import all pipeline modules
// VI: Nhập tất cả module pipeline
// ------------------------------------------------------------
include { VALIDATE_SAMPLESHEET } from './modules/validate_samplesheet'
include { FASTP_QC              } from './modules/qc_fastp'
include { BOWTIE2_ALIGN         } from './modules/align_bowtie2'
include { BAM_FILTER            } from './modules/bam_filter'
include { BIN_COUNT             } from './modules/bin_count'
include { GC_CORRECT            } from './modules/gc_correct'
include { FF_ESTIMATE           } from './modules/ff_estimation'
include { ANEUPLOIDY_CALL       } from './modules/aneuploidy_call'
include { SCA_CALL              } from './modules/sca_call'
include { RAA_SCAN              } from './modules/raa_scan'
include { MICRODELETION_CALL    } from './modules/microdeletion'
include { QC_GATE               } from './modules/qc_gate'
include { MAKE_REPORT           } from './modules/report'

// ------------------------------------------------------------
// EN: Main workflow definition
// VI: Định nghĩa quy trình làm việc chính
// ------------------------------------------------------------
workflow {

    // ── Step 0: Parameter validation / Kiểm tra tham số ─────
    if (!params.samplesheet) {
        error "ERROR / LỖI: --samplesheet is required / --samplesheet là bắt buộc."
    }

    log.info """
    ====================================================
     TTA_NIPTsoft-clinical  v${params.pipeline_version}
     EN: Reference build  : ${params.genome_build}
     EN: Config tier      : ${params.clinical_tier}
     EN: Samplesheet      : ${params.samplesheet}
     VI: Phiên bản hệ gen : ${params.genome_build}
     VI: Cấp độ lâm sàng  : ${params.clinical_tier}
     VI: Samplesheet      : ${params.samplesheet}
    ====================================================
    """.stripIndent()

    // ── Step 1: Validate & parse samplesheet ─────────────────
    // EN: Validate samplesheet CSV against JSON schema; hard-fail on errors
    // VI: Kiểm tra CSV samplesheet dựa trên JSON schema; dừng ngay khi có lỗi
    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)
    VALIDATE_SAMPLESHEET(ch_samplesheet)

    // EN: Parse validated samplesheet into (meta, fastq) channel tuples
    // VI: Phân tích samplesheet đã kiểm tra thành kênh dữ liệu (meta, fastq)
    ch_reads = VALIDATE_SAMPLESHEET.out.validated_csv
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id        : row.sample_id,
                sex       : row.reported_sex ?: 'unknown',
                run_id    : row.run_id ?: params.run_id,
                ga_weeks  : row.ga_weeks ? row.ga_weeks.toFloat() : null,
                batch     : row.batch ?: 'default'
            ]
            def fastq = file(row.fastq_path, checkIfExists: true)
            [ meta, fastq ]
        }

    // ── Step 2: Read QC & trimming / Kiểm tra & cắt tỉa đọc ──
    // EN: fastp filters by quality (Q20), length, adapters, complexity
    // VI: fastp lọc theo chất lượng (Q20), độ dài, adapter, độ phức tạp
    FASTP_QC(ch_reads)

    // ── Step 3: Alignment / Căn chỉnh đọc ────────────────────
    // EN: Bowtie2 single-end alignment, convert and sort to BAM
    // VI: Căn chỉnh Bowtie2 single-end, chuyển đổi và sắp xếp thành BAM
    ch_bowtie2_index = Channel.value(
        file(params.bowtie2_index, checkIfExists: true)
    )
    BOWTIE2_ALIGN(FASTP_QC.out.trimmed_reads, ch_bowtie2_index)

    // ── Step 4: Post-alignment filtering / Lọc sau căn chỉnh ─
    // EN: Remove low-MAPQ, duplicates, blacklist reads
    // VI: Loại bỏ đọc MAPQ thấp, trùng lặp, nằm trong vùng blacklist
    ch_blacklist = Channel.value(
        file(params.blacklist_bed, checkIfExists: true)
    )
    BAM_FILTER(BOWTIE2_ALIGN.out.bam, ch_blacklist)

    // ── Step 5: Genome bin counting / Đếm bins hệ gen ────────
    // EN: Count reads in fixed-size genomic bins (default 50 kb)
    // VI: Đếm số đọc trong các bins hệ gen kích thước cố định (mặc định 50 kb)
    ch_bins_bed = Channel.value(
        file(params.bins_bed, checkIfExists: true)
    )
    BIN_COUNT(BAM_FILTER.out.filtered_bam, ch_bins_bed)

    // ── Step 6: GC / mappability correction / Hiệu chỉnh GC ──
    // EN: LOESS GC correction + mappability normalization → Parquet
    // VI: Hiệu chỉnh LOESS GC + chuẩn hóa tính ánh xạ → file Parquet
    ch_gc_map = Channel.value([
        file(params.gc_content_bed,  checkIfExists: true),
        file(params.mappability_bed, checkIfExists: true)
    ])
    GC_CORRECT(BIN_COUNT.out.raw_counts, ch_gc_map)

    // ── Step 7: Fetal fraction estimation / Ước tính FF ──────
    // EN: SeqFF-style (sex-independent) + Y-based methods
    // VI: Phương pháp kiểu SeqFF (không phụ thuộc giới tính) + dựa trên chrY
    FF_ESTIMATE(GC_CORRECT.out.norm_counts, BAM_FILTER.out.bam_stats)

    // ── Step 8: Aneuploidy calling T21/T18/T13 ───────────────
    // EN: Robust Z-score / NCV model against euploid reference panel
    // VI: Mô hình Z-score bền vững / NCV so với bảng tham chiếu lưỡng bội
    ch_euploid_ref = Channel.value(
        file(params.euploid_reference, checkIfExists: true)
    )
    ANEUPLOIDY_CALL(
        GC_CORRECT.out.norm_counts
            .join(FF_ESTIMATE.out.ff_metrics, by: 0),
        ch_euploid_ref
    )

    // ── Step 9: SCA calling chrX/chrY ────────────────────────
    // EN: Sex chromosome aneuploidy (Turner, Klinefelter, XXX, XYY)
    // VI: Lệch bội NST giới tính (Turner, Klinefelter, XXX, XYY)
    SCA_CALL(
        GC_CORRECT.out.norm_counts
            .join(FF_ESTIMATE.out.ff_metrics, by: 0),
        ch_euploid_ref
    )

    // ── Step 10: RAA scan (tier2 only) / Quét RAA (chỉ tier2) ─
    // EN: Rare autosomal aneuploidy scan with stricter Z ≥ 4.0 threshold
    // VI: Quét lệch bội NST thường hiếm gặp với ngưỡng Z ≥ 4.0 nghiêm hơn
    if (params.clinical_tier == 'tier2' || params.enable_raa) {
        RAA_SCAN(
            GC_CORRECT.out.norm_counts
                .join(FF_ESTIMATE.out.ff_metrics, by: 0),
            ch_euploid_ref
        )
        ch_raa_results = RAA_SCAN.out.raa_calls
    } else {
        ch_raa_results = Channel.empty()
    }

    // ── Step 11: Microdeletion panel (tier2 only) ────────────
    // EN: Panel-constrained microdeletion calling; strictest FF ≥ 8%
    // VI: Gọi vi mất đoạn giới hạn trong panel; FF tối thiểu nghiêm nhất ≥ 8%
    if (params.clinical_tier == 'tier2' || params.enable_microdeletion) {
        ch_microdeletion_panel = Channel.value(
            file(params.microdeletion_panel_bed, checkIfExists: true)
        )
        MICRODELETION_CALL(
            GC_CORRECT.out.norm_counts
                .join(FF_ESTIMATE.out.ff_metrics, by: 0),
            ch_euploid_ref,
            ch_microdeletion_panel
        )
        ch_microdeletion_results = MICRODELETION_CALL.out.microdeletion_calls
    } else {
        ch_microdeletion_results = Channel.empty()
    }

    // ── Step 12: QC gate & no-call decision ──────────────────
    // EN: All gates must pass for a result to be reportable
    // VI: Tất cả cổng QC phải đạt mới được phép báo cáo kết quả
    ch_aneuploidy_calls = ANEUPLOIDY_CALL.out.aneuploidy_calls
    ch_sca_calls        = SCA_CALL.out.sca_calls
    ch_ff_metrics       = FF_ESTIMATE.out.ff_metrics
    ch_fastp_json       = FASTP_QC.out.fastp_json
    ch_bam_stats        = BAM_FILTER.out.bam_stats

    QC_GATE(
        ch_fastp_json
            .join(ch_bam_stats,        by: 0)
            .join(ch_ff_metrics,       by: 0)
            .join(ch_aneuploidy_calls, by: 0)
            .join(ch_sca_calls,        by: 0)
    )

    // ── Step 13: Clinical report generation ──────────────────
    // EN: JSON + HTML report + audit trail for each sample
    // VI: Báo cáo JSON + HTML + dấu kiểm toán cho mỗi mẫu
    ch_qc_decisions = QC_GATE.out.qc_decisions

    MAKE_REPORT(
        ch_qc_decisions
            .join(ch_aneuploidy_calls, by: 0)
            .join(ch_sca_calls,        by: 0)
    )

}

// ------------------------------------------------------------
// EN: Workflow completion handler
// VI: Xử lý sự kiện hoàn thành quy trình làm việc
// ------------------------------------------------------------
workflow.onComplete {
    def status = workflow.success ? "SUCCESS / THÀNH CÔNG" : "FAILED / THẤT BẠI"
    log.info """
    ============================================================
     EN: Pipeline finished  : ${status}
     EN: Completed at       : ${workflow.complete}
     EN: Duration           : ${workflow.duration}
     EN: Work dir           : ${workflow.workDir}
     EN: Exit code          : ${workflow.exitStatus}
     VI: Pipeline kết thúc  : ${status}
     VI: Thời gian hoàn thành: ${workflow.complete}
     VI: Thời gian chạy     : ${workflow.duration}
     VI: Thư mục làm việc   : ${workflow.workDir}
     VI: Mã thoát           : ${workflow.exitStatus}
    ============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "EN: Pipeline error / VI: Lỗi pipeline: ${workflow.errorMessage}"
}
