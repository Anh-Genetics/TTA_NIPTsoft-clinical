// ============================================================
// modules/sca_call.nf
// EN: Sex chromosome aneuploidy (SCA) calling
//     Handles X0 (Turner), XXX (Triple X), XXY (Klinefelter), XYY
// VI: Gọi lệch bội nhiễm sắc thể giới tính (SCA)
//     Xử lý X0 (Turner), XXX, XXY (Klinefelter), XYY
// ============================================================

process SCA_CALL {

    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/calls", mode: 'copy',
        pattern: '*.sca_calls.json'

    conda 'conda-forge::python=3.10 conda-forge::pandas=2.0 conda-forge::scipy=1.11 conda-forge::pyarrow=12.0'

    input:
    tuple val(meta), path(norm_counts), path(ff_metrics)
    path  euploid_reference

    output:
    tuple val(meta), path("${meta.id}.sca_calls.json"), emit: sca_calls

    script:
    def sample_id    = meta.id
    def reported_sex = meta.sex ?: 'unknown'
    def sca_thresh   = params.sca_zscore_threshold ?: 3.0
    def min_ff_xy    = params.sca_min_ff_xy        ?: 0.04

    """
    # EN: Call sex chromosome aneuploidies with sex-specific references
    # VI: Gọi lệch bội NST giới tính với tham chiếu phân tầng theo giới tính
    python3 ${projectDir}/bin/sca_model.py \\
        --norm-counts      ${norm_counts} \\
        --ff-metrics       ${ff_metrics} \\
        --euploid-ref      ${euploid_reference} \\
        --sample-id        ${sample_id} \\
        --reported-sex     ${reported_sex} \\
        --zscore-threshold ${sca_thresh} \\
        --min-ff           ${min_ff_xy} \\
        --out              ${sample_id}.sca_calls.json
    """
}
