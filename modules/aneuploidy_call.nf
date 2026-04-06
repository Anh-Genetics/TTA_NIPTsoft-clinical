// ============================================================
// modules/aneuploidy_call.nf
// Aneuploidy calling for T21, T18, T13
// Uses robust Z-score and NCV models
// ============================================================

process ANEUPLOIDY_CALL {

    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/calls", mode: 'copy',
        pattern: '*.aneuploidy_calls.json'

    conda 'conda-forge::python=3.10 conda-forge::pandas=2.0 conda-forge::scipy=1.11 conda-forge::pyarrow=12.0'

    input:
    tuple val(meta), path(norm_counts), path(ff_metrics)
    path  euploid_reference

    output:
    tuple val(meta), path("${meta.id}.aneuploidy_calls.json"), emit: aneuploidy_calls

    script:
    def sample_id        = meta.id
    def zscore_high      = params.zscore_high_risk   ?: 3.0
    def zscore_grey      = params.zscore_grey_zone_lo ?: 2.5
    def targets          = ['chr21', 'chr18', 'chr13'].join(',')

    """
    python3 ${projectDir}/bin/zscore_model.py \\
        --norm-counts      ${norm_counts} \\
        --ff-metrics       ${ff_metrics} \\
        --euploid-ref      ${euploid_reference} \\
        --sample-id        ${sample_id} \\
        --targets          ${targets} \\
        --zscore-high-risk ${zscore_high} \\
        --zscore-grey-zone ${zscore_grey} \\
        --out              ${sample_id}.aneuploidy_calls.json
    """
}
