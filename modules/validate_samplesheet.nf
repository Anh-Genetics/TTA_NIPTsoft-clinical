// ============================================================
// modules/validate_samplesheet.nf
// Validate CSV samplesheet against JSON schema
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
    python3 ${projectDir}/bin/validate_samplesheet.py \\
        --samplesheet ${samplesheet} \\
        --schema      ${schema} \\
        --out-csv     validated_samplesheet.csv \\
        --out-report  validation_report.json
    """
}
