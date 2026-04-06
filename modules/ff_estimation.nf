// ============================================================
// modules/ff_estimation.nf
// Fetal fraction estimation using sex-independent and optional
// Y-based methods
// ============================================================

process FF_ESTIMATE {

    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/stats", mode: 'copy',
        pattern: '*.ff_metrics.json'

    conda 'conda-forge::python=3.10 conda-forge::pandas=2.0 conda-forge::scikit-learn=1.3 conda-forge::pyarrow=12.0'

    input:
    tuple val(meta), path(norm_counts)
    tuple val(meta2), path(bam_stats)

    output:
    tuple val(meta), path("${meta.id}.ff_metrics.json"), emit: ff_metrics

    script:
    def sample_id     = meta.id
    def reported_sex  = meta.sex ?: 'unknown'

    """
    python3 ${projectDir}/bin/ff_estimate.py \\
        --norm-counts  ${norm_counts} \\
        --bam-stats    ${bam_stats} \\
        --sample-id    ${sample_id} \\
        --reported-sex ${reported_sex} \\
        --out          ${sample_id}.ff_metrics.json
    """
}
