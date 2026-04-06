// ============================================================
// modules/report.nf
// Generate machine-readable JSON clinical report and HTML summary
// ============================================================

process MAKE_REPORT {

    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/reports",  mode: 'copy', pattern: '*.{json,html}'
    publishDir "${params.outdir}/audit",    mode: 'copy', pattern: '*.audit.json'

    conda 'conda-forge::python=3.10 conda-forge::pandas=2.0 conda-forge::jinja2=3.1'

    input:
    tuple val(meta),
          path(qc_decision),
          path(aneuploidy_calls),
          path(sca_calls)

    output:
    tuple val(meta), path("${meta.id}_clinical_report.json"), emit: report_json
    tuple val(meta), path("${meta.id}_clinical_report.html"), emit: report_html
    tuple val(meta), path("${meta.id}.audit.json"),           emit: audit_json

    script:
    def sample_id        = meta.id
    def pipeline_version = params.pipeline_version  ?: '1.0.0'
    def genome_build     = params.genome_build       ?: 'hg19'
    def clinical_tier    = params.clinical_tier      ?: 'tier1'
    def config_version   = params.config_version     ?: '1.0.0'
    def run_id           = meta.run_id               ?: params.run_id

    """
    python3 ${projectDir}/bin/make_report.py \\
        --sample-id        ${sample_id} \\
        --run-id           ${run_id} \\
        --qc-decision      ${qc_decision} \\
        --aneuploidy-calls ${aneuploidy_calls} \\
        --sca-calls        ${sca_calls} \\
        --pipeline-version ${pipeline_version} \\
        --genome-build     ${genome_build} \\
        --clinical-tier    ${clinical_tier} \\
        --config-version   ${config_version} \\
        --report-schema    ${projectDir}/schemas/report.schema.json \\
        --out-json         ${sample_id}_clinical_report.json \\
        --out-html         ${sample_id}_clinical_report.html \\
        --out-audit        ${sample_id}.audit.json
    """
}
