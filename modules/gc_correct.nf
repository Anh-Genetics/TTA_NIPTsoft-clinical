// ============================================================
// modules/gc_correct.nf
// GC/mappability-aware normalization using LOESS regression
// ============================================================

process GC_CORRECT {

    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/counts", mode: 'copy',
        pattern: '*.norm_counts.parquet'
    publishDir "${params.outdir}/qc/gc",  mode: 'copy',
        pattern: '*.gc_bias.json'

    conda 'conda-forge::r-base=4.3 bioconda::bioconductor-genomicranges=1.54 conda-forge::r-arrow=12.0'

    input:
    tuple val(meta), path(raw_counts)
    tuple path(gc_bed), path(mappability_bed)

    output:
    tuple val(meta), path("${meta.id}.norm_counts.parquet"), emit: norm_counts
    tuple val(meta), path("${meta.id}.gc_bias.json"),         emit: gc_bias_metrics

    script:
    def sample_id = meta.id

    """
    Rscript ${projectDir}/bin/gc_loess_correct.R \\
        --counts       ${raw_counts} \\
        --gc-bed       ${gc_bed} \\
        --map-bed      ${mappability_bed} \\
        --sample-id    ${sample_id} \\
        --out-parquet  ${sample_id}.norm_counts.parquet \\
        --out-metrics  ${sample_id}.gc_bias.json
    """
}
