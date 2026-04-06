// ============================================================
// modules/bin_count.nf
// EN: Count reads per genomic bin using bedtools coverage
// VI: Đếm số đọc trong mỗi bin hệ gen bằng bedtools coverage
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
    echo "[INFO] EN: Counting reads per bin / VI: Đang đếm số đọc trong mỗi bin cho: ${sample_id}" >&2

    # EN: Count reads overlapping each bin in the BED file
    # VI: Đếm số đọc chồng lên mỗi bin trong file BED
    # EN: Output columns: chr start end name gc_content mappability raw_count
    # VI: Các cột đầu ra: chr start end tên gc_content tính_ánh_xạ số_đọc
    bedtools coverage \\
        -a    ${bins_bed} \\
        -b    ${bam} \\
        -counts \\
    > ${sample_id}.counts.bed

    # EN: Verify output is non-empty — fail fast if no bins were counted
    # VI: Kiểm tra đầu ra không trống — dừng ngay nếu không có bin nào được đếm
    if [ ! -s ${sample_id}.counts.bed ]; then
        echo "ERROR / LỖI: Empty bin counts / Số đếm bins trống cho mẫu: ${sample_id}" >&2
        exit 1
    fi

    echo "[DONE / HOÀN THÀNH] EN: Bin counting complete / VI: Hoàn thành đếm bins cho: ${sample_id}" >&2
    """
}
