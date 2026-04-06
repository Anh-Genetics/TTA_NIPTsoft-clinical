// ============================================================
// modules/bam_filter.nf
// Post-alignment filtering:
//   - Remove blacklist / low-mappability regions
//   - Apply MAPQ threshold
//   - Mark duplicates (configurable remove)
//   - Collect alignment metrics
// ============================================================

process BAM_FILTER {

    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/bam",       mode: 'copy',
        pattern: '*.filtered.bam*',  enabled: params.publish_bam
    publishDir "${params.outdir}/qc/bam",    mode: 'copy',
        pattern: '*.stats.txt'

    conda 'bioconda::samtools=1.17 bioconda::picard=3.1.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path  blacklist_bed

    output:
    tuple val(meta), path("${meta.id}.filtered.bam"),     emit: filtered_bam
    tuple val(meta), path("${meta.id}.filtered.bam.bai"), emit: filtered_bai
    tuple val(meta), path("${meta.id}.bam_stats.txt"),    emit: bam_stats
    tuple val(meta), path("${meta.id}.dup_metrics.txt"),  emit: dup_metrics

    script:
    def sample_id      = meta.id
    def mapq           = params.mapq_threshold    ?: 30
    def mark_dup       = params.mark_duplicates   != false
    def remove_dup     = params.remove_duplicates ? 'true' : 'false'
    def dup_flag       = params.remove_duplicates ? '-F 1024' : ''
    def threads        = task.cpus

    """
    # Step 1: Remove blacklisted regions and apply MAPQ filter
    samtools view -b -q ${mapq} -F 4 -F 256 -F 2048 ${bam} \\
    | bedtools intersect -v -abam stdin -b ${blacklist_bed} \\
    | samtools sort -@ ${threads} -o tmp_mapq_filtered.bam -

    samtools index tmp_mapq_filtered.bam

    # Step 2: Mark (and optionally remove) duplicates with Picard
    if [ "${mark_dup}" = "true" ]; then
        picard MarkDuplicates \\
            INPUT=tmp_mapq_filtered.bam \\
            OUTPUT=tmp_dedup.bam \\
            METRICS_FILE=${sample_id}.dup_metrics.txt \\
            REMOVE_DUPLICATES=${remove_dup} \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=./tmp_picard \\
            QUIET=true
    else
        cp tmp_mapq_filtered.bam tmp_dedup.bam
        touch ${sample_id}.dup_metrics.txt
    fi

    # Step 3: Sort and index final BAM
    samtools sort -@ ${threads} -o ${sample_id}.filtered.bam tmp_dedup.bam
    samtools index ${sample_id}.filtered.bam

    # Step 4: Generate alignment statistics
    samtools flagstat ${sample_id}.filtered.bam > ${sample_id}.bam_stats.txt
    samtools idxstats  ${sample_id}.filtered.bam >> ${sample_id}.bam_stats.txt

    # Cleanup
    rm -f tmp_mapq_filtered.bam tmp_mapq_filtered.bam.bai tmp_dedup.bam
    """
}
