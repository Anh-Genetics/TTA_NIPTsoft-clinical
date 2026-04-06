// ============================================================
// modules/bin_count.nf
// Count reads per genomic bin using bedtools coverage
// ============================================================

process BIN_COUNT {

    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/counts", mode: 'copy', pattern: '*.counts.bed'

    conda 'bioconda::bedtools=2.31.0 bioconda::samtools=1.17'

    input:
    tuple val(meta), path(bam), path(bai)
    path  bins_bed

    output:
    tuple val(meta), path("${meta.id}.counts.bed"), emit: raw_counts

    script:
    def sample_id = meta.id

    """
    # Count reads overlapping each bin
    # Output: chr start end name gc_content mappability raw_count
    bedtools coverage \\
        -a    ${bins_bed} \\
        -b    ${bam} \\
        -counts \\
    > ${sample_id}.counts.bed

    # Verify output is non-empty
    if [ ! -s ${sample_id}.counts.bed ]; then
        echo "ERROR: Empty bin counts for sample ${sample_id}" >&2
        exit 1
    fi
    """
}
