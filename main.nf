#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nfdata-omics/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nfdata-omics/rnaseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ                  } from './workflows/rnaseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PREPARE_REFERENCES      } from './subworkflows/local/prepare_references'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Use getGenomeAttribute() to fetch parameters from igenomes.config using `--genome`
params.fasta            = getGenomeAttribute('fasta')
params.gtf              = getGenomeAttribute('gtf')
params.transcript_fasta = getGenomeAttribute('transcript_fasta')
params.star_index       = getGenomeAttribute('star')
params.salmon_index     = getGenomeAttribute('salmon')
params.ref_flat         = getGenomeAttribute('ref_flat')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFDATAOMICS_RNASEQ {

    take:
    samplesheet // channel: samplesheet read in from --input
    counts      // channel: count matrix file read as --counts
    metadata    // channel: sample metadata table read as --metadata

    main:

    ch_versions = channel.empty()

    //
    // SUBWORKFLOW: Prepare reference files
    //
    PREPARE_REFERENCES(
        params.fasta,
        params.gtf,
        params.transcript_fasta,
        params.star_index,
        params.salmon_index,
        params.ref_flat
    )
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

    //
    // WORKFLOW: Run pipeline
    //
    RNASEQ (
        samplesheet,
        counts,
        metadata,
        PREPARE_REFERENCES.out.fasta,
        PREPARE_REFERENCES.out.gtf,
        PREPARE_REFERENCES.out.fai,
        PREPARE_REFERENCES.out.chrom_sizes,
        PREPARE_REFERENCES.out.gene_bed,
        PREPARE_REFERENCES.out.transcript_fasta,
        PREPARE_REFERENCES.out.star_index,
        PREPARE_REFERENCES.out.salmon_index,
        PREPARE_REFERENCES.out.refflat,
        ch_versions,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir
    )
    ch_versions = ch_versions.mix(RNASEQ.out.versions)

    emit:
    multiqc_report = RNASEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions               // channel: [version1, version2, ...]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.counts,
        params.metadata,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFDATAOMICS_RNASEQ (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.counts,
        PIPELINE_INITIALISATION.out.metadata
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFDATAOMICS_RNASEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
