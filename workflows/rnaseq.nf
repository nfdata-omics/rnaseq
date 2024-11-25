/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 }  from '../modules/nf-core/fastqc/main'
include { MULTIQC                }  from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       }  from 'plugin/nf-schema'
include { paramsSummaryMultiqc   }  from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML }  from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText }  from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { DIFFERENTIAL_EXPRESSION } from '../subworkflows/local/differential_expression'
include { DIM_REDUCTION }           from '../subworkflows/local/dim_reduction'

include { PICARD_COLLECTRNASEQMETRICS } from '../modules/nf-core/picard/collectrnaseqmetrics/main'
include { SAMPLE_FILTER }               from '../modules/local/sample_filter'
include { SPLIT_COUNT_MATRIX }          from '../modules/local/split_count_matrix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    counts         // channel: count matrix file read as --counts
    metadata       // channel: sample metadata table read as --metadata

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet
        .filter{ it[0].data_type == "fastq" }
        .map {
            meta, fastq1s, fastq2s, _bams ->
                return [ meta, fastq1s + fastq2s ]
        }
        .set { ch_fastqs }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastqs
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    ch_samplesheet
        .filter{ it[0].data_type == "bam" }
        .map {
            meta, _fastq1s, _fastq2s, bams ->
                return [ meta, bams ]
        }
        .set { ch_bam }

    //
    // INSERT BAM PROCESSING HERE
    //

    PICARD_COLLECTRNASEQMETRICS(
        ch_bam,
        params.ref_flat,
        params.fasta,
        []
    )

    //
    // DOWNSTREAM ANALYSIS OF COUNT MATRIX
    //

    if ( params.exclude_list ) {
        SAMPLE_FILTER(
            counts,
            params.exclude_list
        )
        counts_filt = SAMPLE_FILTER.out.filtered
        ch_versions = ch_versions.mix(SAMPLE_FILTER.out.versions)
    } else {
        counts_filt = counts
    }

    //
    // DIMENSIONALITY REDUCTION SUBWORKFLOW
    //

    counts_filt
        .map { file -> [["id":file.baseName], file] }
        .set { ch_counts }

    DIM_REDUCTION(
        ch_counts,
        params.gene_column_nr,
        params.gene_id_index,
        metadata,
        params.frac_expressed
    )
    ch_versions = ch_versions.mix(DIM_REDUCTION.out.versions)

    // COUNT SUBSETS

    if ( params.split_file ) {
        SPLIT_COUNT_MATRIX(
            counts_filt,
            params.gene_column_nr,
            params.split_file
        )
        counts_split = SPLIT_COUNT_MATRIX.out.matrices.flatten()
        ch_versions = ch_versions.mix(SPLIT_COUNT_MATRIX.out.versions)
    } else {
        counts_split = counts_filt
    }

    counts_split
        .map { file -> [["id":file.baseName], file] }
        .set { ch_counts }

    // DEA AND FUNCTIONAL

    DIFFERENTIAL_EXPRESSION(
        ch_counts,
        params.gene_column_nr,
        params.gene_id_index,
        metadata,
        params.model_formula,
        params.comparisons,
        params.frac_expressed,
        params.fdr_threshold,
        params.lfc_threshold,
    )
    ch_versions = ch_versions.mix(DIFFERENTIAL_EXPRESSION.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
