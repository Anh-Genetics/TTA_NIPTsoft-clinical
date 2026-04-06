// ============================================================
// modules/align_bowtie2.nf
// Align SE reads with Bowtie2, convert to sorted BAM
// ============================================================

process BOWTIE2_ALIGN {

    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/bam", mode: 'copy',
        pattern: '*.bam*', enabled: params.publish_bam

    conda 'bioconda::bowtie2=2.5.1 bioconda::samtools=1.17'

    input:
    tuple val(meta), path(reads)
    path  bowtie2_index_dir

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.sorted.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}_align.log"),      emit: align_log

    script:
    def sample_id = meta.id
    def idx_base  = "${bowtie2_index_dir}/${bowtie2_index_dir.getName()}"
    // Prefer bowtie2_index param-based path when index is a directory
    def bt2_idx   = params.bowtie2_index
    def mode      = params.bowtie2_mode    ?: '--sensitive'
    def no_unal   = params.bowtie2_no_unal ? '--no-unal' : ''
    def threads   = task.cpus

    """
    bowtie2 \\
        -x  ${bt2_idx} \\
        -U  ${reads} \\
        ${mode} \\
        ${no_unal} \\
        --threads ${threads} \\
        2>${sample_id}_align.log \\
    | samtools view -bS -q 1 - \\
    | samtools sort -@ ${threads} -o ${sample_id}.sorted.bam -

    samtools index ${sample_id}.sorted.bam
    """
}
