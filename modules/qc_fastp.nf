// ============================================================
// modules/qc_fastp.nf
// Read QC, trimming, and filtering with fastp (SE mode)
// ============================================================

process FASTP_QC {

    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/qc/fastp",  mode: 'copy', pattern: '*.{json,html}'
    publishDir "${params.outdir}/reads",     mode: 'copy', pattern: '*_trimmed.fastq.gz', enabled: false

    conda 'bioconda::fastp=0.23.4'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(meta), path("${meta.id}_fastp.json"),        emit: fastp_json
    tuple val(meta), path("${meta.id}_fastp.html"),        emit: fastp_html

    script:
    def sample_id         = meta.id
    def min_length        = params.min_read_length        ?: 30
    def quality_threshold = params.quality_threshold      ?: 20
    def adapter           = params.adapter_sequence       ?: ''
    def adapter_arg       = adapter ? "--adapter_sequence ${adapter}" : '--detect_adapter_for_pe'
    def complexity_filter = (params.fastp_low_complexity_filter) ? "--low_complexity_filter --complexity_threshold ${params.fastp_complexity_threshold ?: 30}" : ''
    def n_base_limit      = params.fastp_n_base_limit     ?: 5
    def avg_qual          = params.fastp_average_qual     ?: 20
    def threads           = task.cpus

    """
    fastp \\
        --in1           ${reads} \\
        --out1          ${sample_id}_trimmed.fastq.gz \\
        ${adapter_arg} \\
        --length_required   ${min_length} \\
        --average_qual      ${avg_qual} \\
        --n_base_limit      ${n_base_limit} \\
        ${complexity_filter} \\
        --disable_quality_filtering false \\
        --qualified_quality_phred ${quality_threshold} \\
        --thread            ${threads} \\
        --json              ${sample_id}_fastp.json \\
        --html              ${sample_id}_fastp.html \\
        2>&1 | tee ${sample_id}_fastp.log
    """
}
