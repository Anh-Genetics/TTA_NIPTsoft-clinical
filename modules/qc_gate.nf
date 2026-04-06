// ============================================================
// modules/qc_gate.nf
// EN: Unified QC gate and no-call decision engine.
//     Aggregates all QC metrics; determines final reportability.
//     All gates (G01–G06) must pass for results to be reported.
// VI: Bộ máy kiểm tra QC thống nhất và quyết định không gọi kết quả.
//     Tổng hợp tất cả chỉ số QC; xác định khả năng báo cáo cuối cùng.
//     Tất cả cổng (G01–G06) phải đạt mới được phép báo cáo kết quả.
// ============================================================

process QC_GATE {

    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/qc", mode: 'copy',
        pattern: '*.qc_decision.json'

    conda 'conda-forge::python=3.10 conda-forge::pandas=2.0'

    input:
    tuple val(meta),
          path(fastp_json),
          path(bam_stats),
          path(ff_metrics),
          path(aneuploidy_calls),
          path(sca_calls)

    output:
    tuple val(meta), path("${meta.id}.qc_decision.json"), emit: qc_decisions

    script:
    def sample_id      = meta.id
    // EN: QC thresholds with safe defaults matching clinical_tier1.config
    // VI: Ngưỡng QC với giá trị mặc định an toàn khớp với clinical_tier1.config
    def min_raw_reads  = params.min_raw_reads      ?: 8000000
    def min_uniq_reads = params.min_unique_reads   ?: 4000000
    def max_dup_rate   = params.max_dup_rate        ?: 0.35
    def min_ff         = params.min_fetal_fraction  ?: 0.04
    def min_map_rate   = params.min_mapping_rate    ?: 0.70

    """
    # EN: Evaluate all QC gates and determine reportability
    # VI: Đánh giá tất cả cổng QC và xác định khả năng báo cáo
    python3 ${projectDir}/bin/qc_decision.py \\
        --sample-id        ${sample_id} \\
        --fastp-json       ${fastp_json} \\
        --bam-stats        ${bam_stats} \\
        --ff-metrics       ${ff_metrics} \\
        --aneuploidy-calls ${aneuploidy_calls} \\
        --sca-calls        ${sca_calls} \\
        --min-raw-reads    ${min_raw_reads} \\
        --min-unique-reads ${min_uniq_reads} \\
        --max-dup-rate     ${max_dup_rate} \\
        --min-ff           ${min_ff} \\
        --min-mapping-rate ${min_map_rate} \\
        --out              ${sample_id}.qc_decision.json
    """
}
