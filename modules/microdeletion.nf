// ============================================================
// modules/microdeletion.nf
// EN: Microdeletion detection constrained to curated panel regions.
//     Strictest QC requirements (FF ≥ 8%); conservative no-call behavior.
//     HIGH_RISK result MUST be confirmed by diagnostic testing.
// VI: Phát hiện vi mất đoạn giới hạn trong panel đã tuyển chọn.
//     Yêu cầu QC nghiêm nhất (FF ≥ 8%); hành vi thận trọng khi không gọi.
//     Kết quả CAO NGUY CƠ PHẢI được xác nhận bằng xét nghiệm chẩn đoán.
// ============================================================

process MICRODELETION_CALL {

    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/calls", mode: 'copy',
        pattern: '*.microdeletion_calls.json'

    conda 'conda-forge::python=3.10 conda-forge::pandas=2.0 conda-forge::scipy=1.11 conda-forge::pyarrow=12.0'

    input:
    tuple val(meta), path(norm_counts), path(ff_metrics)
    path  euploid_reference
    path  panel_bed

    output:
    tuple val(meta), path("${meta.id}.microdeletion_calls.json"), emit: microdeletion_calls

    script:
    def sample_id   = meta.id
    def md_thresh   = params.microdeletion_zscore   ?: 5.0
    def md_min_ff   = params.microdeletion_min_ff   ?: 0.08
    def md_min_bins = params.microdeletion_min_bins ?: 2

    """
    # EN: Detect panel-region microdeletions (deletion = negative Z ≤ -5)
    # VI: Phát hiện vi mất đoạn trong vùng panel (mất đoạn = Z âm ≤ -5)
    python3 ${projectDir}/bin/microdeletion_call.py \\
        --norm-counts      ${norm_counts} \\
        --ff-metrics       ${ff_metrics} \\
        --euploid-ref      ${euploid_reference} \\
        --panel-bed        ${panel_bed} \\
        --sample-id        ${sample_id} \\
        --zscore-threshold ${md_thresh} \\
        --min-ff           ${md_min_ff} \\
        --min-bins         ${md_min_bins} \\
        --out              ${sample_id}.microdeletion_calls.json
    """
}
