#!/usr/bin/env nextflow
// ============================================================
// TTA_NIPTsoft-clinical  — Illumina single-end NIPT pipeline
// Nextflow DSL2 main orchestrator
// ============================================================
nextflow.enable.dsl = 2

// ------------------------------------------------------------
// Import modules
// ------------------------------------------------------------
include { VALIDATE_SAMPLESHEET } from './modules/validate_samplesheet'
include { FASTP_QC              } from './modules/qc_fastp'
include { BOWTIE2_ALIGN         } from './modules/align_bowtie2'
include { BAM_FILTER            } from './modules/bam_filter'
include { BIN_COUNT             } from './modules/bin_count'
include { GC_CORRECT            } from './modules/gc_correct'
include { FF_ESTIMATE           } from './modules/ff_estimation'
include { ANEUPLOIDY_CALL       } from './modules/aneuploidy_call'
include { SCA_CALL              } from './modules/sca_call'
include { RAA_SCAN              } from './modules/raa_scan'
include { MICRODELETION_CALL    } from './modules/microdeletion'
include { QC_GATE               } from './modules/qc_gate'
include { MAKE_REPORT           } from './modules/report'

// ------------------------------------------------------------
// Pipeline workflow
// ------------------------------------------------------------
workflow {

    // ── 0. Parameter validation ──────────────────────────────
    if (!params.samplesheet) {
        error "ERROR: --samplesheet is required."
    }

    log.info """
    ====================================================
     TTA_NIPTsoft-clinical  v${params.pipeline_version}
     Reference build : ${params.genome_build}
     Config tier     : ${params.clinical_tier}
     Samplesheet     : ${params.samplesheet}
    ====================================================
    """.stripIndent()

    // ── 1. Validate and parse samplesheet ───────────────────
    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)

    VALIDATE_SAMPLESHEET(ch_samplesheet)

    // Parse validated samplesheet into (meta, fastq) tuples
    ch_reads = VALIDATE_SAMPLESHEET.out.validated_csv
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id        : row.sample_id,
                sex       : row.reported_sex ?: 'unknown',
                run_id    : row.run_id ?: params.run_id,
                ga_weeks  : row.ga_weeks ? row.ga_weeks.toFloat() : null,
                batch     : row.batch ?: 'default'
            ]
            def fastq = file(row.fastq_path, checkIfExists: true)
            [ meta, fastq ]
        }

    // ── 2. Read QC & trimming ────────────────────────────────
    FASTP_QC(ch_reads)

    // ── 3. Alignment ─────────────────────────────────────────
    ch_bowtie2_index = Channel.value(
        file(params.bowtie2_index, checkIfExists: true)
    )
    BOWTIE2_ALIGN(FASTP_QC.out.trimmed_reads, ch_bowtie2_index)

    // ── 4. Post-alignment filtering ──────────────────────────
    ch_blacklist = Channel.value(
        file(params.blacklist_bed, checkIfExists: true)
    )
    BAM_FILTER(BOWTIE2_ALIGN.out.bam, ch_blacklist)

    // ── 5. Bin counting ──────────────────────────────────────
    ch_bins_bed = Channel.value(
        file(params.bins_bed, checkIfExists: true)
    )
    BIN_COUNT(BAM_FILTER.out.filtered_bam, ch_bins_bed)

    // ── 6. GC / mappability correction ───────────────────────
    ch_gc_map = Channel.value([
        file(params.gc_content_bed,      checkIfExists: true),
        file(params.mappability_bed,     checkIfExists: true)
    ])
    GC_CORRECT(BIN_COUNT.out.raw_counts, ch_gc_map)

    // ── 7. Fetal fraction estimation ─────────────────────────
    FF_ESTIMATE(GC_CORRECT.out.norm_counts, BAM_FILTER.out.bam_stats)

    // ── 8. Aneuploidy calling (T21 / T18 / T13) ──────────────
    ch_euploid_ref = Channel.value(
        file(params.euploid_reference, checkIfExists: true)
    )
    ANEUPLOIDY_CALL(
        GC_CORRECT.out.norm_counts
            .join(FF_ESTIMATE.out.ff_metrics, by: 0),
        ch_euploid_ref
    )

    // ── 9. SCA calling ───────────────────────────────────────
    SCA_CALL(
        GC_CORRECT.out.norm_counts
            .join(FF_ESTIMATE.out.ff_metrics, by: 0),
        ch_euploid_ref
    )

    // ── 10. RAA scan (if tier2) ──────────────────────────────
    if (params.clinical_tier == 'tier2' || params.enable_raa) {
        RAA_SCAN(
            GC_CORRECT.out.norm_counts
                .join(FF_ESTIMATE.out.ff_metrics, by: 0),
            ch_euploid_ref
        )
        ch_raa_results = RAA_SCAN.out.raa_calls
    } else {
        ch_raa_results = Channel.empty()
    }

    // ── 11. Microdeletion panel (if tier2) ───────────────────
    if (params.clinical_tier == 'tier2' || params.enable_microdeletion) {
        ch_microdeletion_panel = Channel.value(
            file(params.microdeletion_panel_bed, checkIfExists: true)
        )
        MICRODELETION_CALL(
            GC_CORRECT.out.norm_counts
                .join(FF_ESTIMATE.out.ff_metrics, by: 0),
            ch_euploid_ref,
            ch_microdeletion_panel
        )
        ch_microdeletion_results = MICRODELETION_CALL.out.microdeletion_calls
    } else {
        ch_microdeletion_results = Channel.empty()
    }

    // ── 12. QC gate & no-call decision ───────────────────────
    ch_aneuploidy_calls = ANEUPLOIDY_CALL.out.aneuploidy_calls
    ch_sca_calls        = SCA_CALL.out.sca_calls
    ch_ff_metrics       = FF_ESTIMATE.out.ff_metrics
    ch_fastp_json       = FASTP_QC.out.fastp_json
    ch_bam_stats        = BAM_FILTER.out.bam_stats

    QC_GATE(
        ch_fastp_json
            .join(ch_bam_stats,        by: 0)
            .join(ch_ff_metrics,       by: 0)
            .join(ch_aneuploidy_calls, by: 0)
            .join(ch_sca_calls,        by: 0)
    )

    // ── 13. Clinical report ───────────────────────────────────
    ch_qc_decisions = QC_GATE.out.qc_decisions

    MAKE_REPORT(
        ch_qc_decisions
            .join(ch_aneuploidy_calls, by: 0)
            .join(ch_sca_calls,        by: 0)
    )

}

// ------------------------------------------------------------
// Workflow completion handler
// ------------------------------------------------------------
workflow.onComplete {
    def status = workflow.success ? "SUCCESS" : "FAILED"
    log.info """
    ============================================================
     Pipeline finished  : ${status}
     Completed at       : ${workflow.complete}
     Duration           : ${workflow.duration}
     Work dir           : ${workflow.workDir}
     Exit code          : ${workflow.exitStatus}
    ============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
}
