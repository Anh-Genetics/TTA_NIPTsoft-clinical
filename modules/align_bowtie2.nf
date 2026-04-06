// ============================================================
// modules/align_bowtie2.nf
// EN: Align single-end reads with Bowtie2 and convert to sorted BAM
// VI: Căn chỉnh đọc single-end bằng Bowtie2 và chuyển đổi thành BAM đã sắp xếp
// ============================================================

process BOWTIE2_ALIGN {

    tag "${meta.id}"
    label 'process_high'

    // EN: Publish BAM files only when params.publish_bam is true (saves disk)
    // VI: Xuất file BAM chỉ khi params.publish_bam là true (tiết kiệm dung lượng)
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
    // EN: Use bowtie2_index param directly (full path to index prefix)
    // VI: Dùng tham số bowtie2_index trực tiếp (đường dẫn đầy đủ đến tiền tố index)
    def bt2_idx   = params.bowtie2_index
    def mode      = params.bowtie2_mode    ?: '--sensitive'
    def no_unal   = params.bowtie2_no_unal ? '--no-unal' : ''
    def threads   = task.cpus

    """
    # EN: Align SE reads → SAM → filter primary alignments → sort → BAM
    # VI: Căn chỉnh đọc SE → SAM → lọc căn chỉnh chính → sắp xếp → BAM
    echo "[INFO] EN: Aligning reads / VI: Đang căn chỉnh đọc cho mẫu: ${sample_id}" >&2

    bowtie2 \\
        -x  ${bt2_idx} \\
        -U  ${reads} \\
        ${mode} \\
        ${no_unal} \\
        --threads ${threads} \\
        2>${sample_id}_align.log \\
    | samtools view -bS -q 1 - \\
    | samtools sort -@ ${threads} -o ${sample_id}.sorted.bam -

    # EN: Index the sorted BAM for downstream tools
    # VI: Tạo chỉ mục cho BAM đã sắp xếp cho các công cụ tiếp theo
    samtools index ${sample_id}.sorted.bam

    echo "[DONE / HOÀN THÀNH] EN: Alignment complete / VI: Hoàn thành căn chỉnh cho: ${sample_id}" >&2
    """
}
