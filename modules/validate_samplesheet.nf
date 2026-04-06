// ============================================================
// modules/validate_samplesheet.nf
// EN: Validate CSV samplesheet against JSON schema.
//     Fails hard on any error (duplicate IDs, missing fields, etc.)
// VI: Kiểm tra file CSV samplesheet dựa trên JSON schema.
//     Dừng ngay khi có lỗi (ID trùng, trường thiếu, v.v.)
// ============================================================

process VALIDATE_SAMPLESHEET {

    tag "samplesheet"
    label 'process_single'

    publishDir "${params.outdir}/audit", mode: 'copy'

    conda 'conda-forge::python=3.10 conda-forge::jsonschema=4.17'

    input:
    path samplesheet

    output:
    path 'validated_samplesheet.csv', emit: validated_csv
    path 'validation_report.json',    emit: validation_report

    script:
    def schema = "${projectDir}/schemas/samplesheet.schema.json"
    """
    # EN: Validate samplesheet CSV against JSON schema; produce cleaned CSV
    # VI: Kiểm tra CSV samplesheet dựa trên JSON schema; tạo CSV đã làm sạch
    python3 ${projectDir}/bin/validate_samplesheet.py \\
        --samplesheet ${samplesheet} \\
        --schema      ${schema} \\
        --out-csv     validated_samplesheet.csv \\
        --out-report  validation_report.json
    """
}
