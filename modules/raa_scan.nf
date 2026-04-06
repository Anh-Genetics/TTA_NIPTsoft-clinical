// ============================================================
// modules/raa_scan.nf
// EN: Rare Autosomal Aneuploidy (RAA) scan across all autosomes
//     Uses stricter Z ≥ 4.0 threshold and FF ≥ 6% requirement.
//     Results are "suspected" flags, not definitive calls.
// VI: Quét lệch bội NST thường hiếm gặp (RAA) trên tất cả NST thường
//     Dùng ngưỡng Z ≥ 4.0 nghiêm hơn và yêu cầu FF ≥ 6%.
//     Kết quả là cờ "nghi ngờ", không phải kết luận dứt khoát.
// ============================================================

process RAA_SCAN {

    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/calls", mode: 'copy',
        pattern: '*.raa_calls.json'

    conda 'conda-forge::python=3.10 conda-forge::pandas=2.0 conda-forge::scipy=1.11 conda-forge::pyarrow=12.0'

    input:
    tuple val(meta), path(norm_counts), path(ff_metrics)
    path  euploid_reference

    output:
    tuple val(meta), path("${meta.id}.raa_calls.json"), emit: raa_calls

    script:
    def sample_id   = meta.id
    def raa_thresh  = params.raa_zscore_threshold ?: 4.0
    def min_bins    = params.raa_min_size_bins    ?: 5
    def min_ff      = params.raa_min_ff           ?: 0.06
    def chromosomes = (params.raa_chromosomes     ?: []).join(',')
    if (!chromosomes) {
        // EN: Default: all normalizer autosomes (T21/T18/T13 included for completeness)
        // VI: Mặc định: tất cả NST thường chuẩn hóa (T21/T18/T13 bao gồm đầy đủ)
        chromosomes = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr14,chr15,chr16,chr17,chr19,chr20,chr22'
    }

    """
    # EN: Scan all specified autosomes for rare aneuploidies
    # VI: Quét tất cả NST thường được chỉ định tìm lệch bội hiếm gặp
    python3 ${projectDir}/bin/raa_scan.py \\
        --norm-counts      ${norm_counts} \\
        --ff-metrics       ${ff_metrics} \\
        --euploid-ref      ${euploid_reference} \\
        --sample-id        ${sample_id} \\
        --chromosomes      ${chromosomes} \\
        --zscore-threshold ${raa_thresh} \\
        --min-bins         ${min_bins} \\
        --min-ff           ${min_ff} \\
        --out              ${sample_id}.raa_calls.json
    """
}
