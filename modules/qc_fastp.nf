// ============================================================
// modules/qc_fastp.nf
// EN: Read QC, adapter trimming, and quality filtering with fastp (SE mode)
// VI: Kiểm tra chất lượng đọc, cắt adapter và lọc chất lượng bằng fastp (chế độ SE)
// ============================================================

process FASTP_QC {

    tag "${meta.id}"
    label 'process_medium'

    // EN: Publish QC reports to outdir/qc/fastp; trimmed reads not kept by default
    // VI: Xuất báo cáo QC vào outdir/qc/fastp; đọc đã cắt tỉa không lưu mặc định
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
    // EN: Quality and length thresholds from params with safe defaults
    // VI: Ngưỡng chất lượng và độ dài từ tham số với giá trị mặc định an toàn
    def min_length        = params.min_read_length        ?: 30
    def quality_threshold = params.quality_threshold      ?: 20
    def adapter           = params.adapter_sequence       ?: ''
    def adapter_arg       = adapter ? "--adapter_sequence ${adapter}" : '--detect_adapter_for_pe'
    def complexity_filter = (params.fastp_low_complexity_filter) ? "--low_complexity_filter --complexity_threshold ${params.fastp_complexity_threshold ?: 30}" : ''
    def n_base_limit      = params.fastp_n_base_limit     ?: 5
    def avg_qual          = params.fastp_average_qual     ?: 20
    def threads           = task.cpus

    """
    # EN: Run fastp for SE read QC, trimming, and quality filtering
    # VI: Chạy fastp để kiểm tra QC, cắt tỉa và lọc chất lượng đọc SE
    echo "[INFO] EN: Starting fastp QC / VI: Bắt đầu kiểm tra fastp cho mẫu: ${sample_id}" >&2

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

    echo "[DONE / HOÀN THÀNH] EN: fastp QC complete / VI: Hoàn thành kiểm tra fastp cho: ${sample_id}" >&2
    """
}
