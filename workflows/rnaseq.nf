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
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS } from '../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
include { DIFFERENTIAL_EXPRESSION } from '../subworkflows/local/differential_expression'
include { DIM_REDUCTION }           from '../subworkflows/local/dim_reduction'
include { BAM_MARKDUPLICATES_PICARD } from '../subworkflows/nf-core/bam_markduplicates_picard/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENES } from '../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_BIOTYPES } from '../modules/nf-core/subread/featurecounts/main'
include { multiqcTsvFromList     } from '../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'

include { PICARD_COLLECTRNASEQMETRICS } from '../modules/nf-core/picard/collectrnaseqmetrics/main'
include { SAMPLE_FILTER }               from '../modules/local/sample_filter'
include { SPLIT_COUNT_MATRIX }          from '../modules/local/split_count_matrix'
include { STAR_ALIGN                  } from '../modules/nf-core/star/align'
include { MULTIQC_CUSTOM_BIOTYPE             } from '../modules/local/multiqc_custom_biotype'
include { GENEID_TO_GENENAME             } from '../modules/local/geneid_to_genename'
include { JOIN_FEATURECOUNTS_MATRIX      } from '../modules/local/join_featurecounts_matrix/main'
include { SAMTOOLS_SORT      } from '../modules/nf-core/samtools/sort/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow RNASEQ {

    take:
    ch_samplesheet       // channel: samplesheet read in from --input
    counts               // channel: count matrix file read as --counts
    metadata             // channel: sample metadata table read as --metadata
    ch_fasta             // channel: path(genome.fasta)
    ch_gtf               // channel: path(genome.gtf)
    ch_fai               // channel: path(genome.fai)
    ch_chrom_sizes       // channel: path(genome.sizes)
    ch_gene_bed          // channel: path(gene.bed)
    ch_transcript_fasta  // channel: path(transcript.fasta)
    ch_star_index        // channel: path(star/index/)
    ch_salmon_index      // channel: path(salmon/index/)
    ch_ref_flat          // channel: path(refFlat.txt)
    ch_versions          // channel: [ path(versions.yml) ]

    main:

    ch_multiqc_files = Channel.empty()

    // Header files for MultiQC
    ch_pca_header_multiqc           = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
    sample_status_header_multiqc    = file("$projectDir/assets/multiqc/sample_status_header.txt", checkIfExists: true)
    ch_clustering_header_multiqc    = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
    ch_biotypes_header_multiqc      = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

    //
    // Create channel from input file provided through params.input
    //
    ch_samplesheet
        .filter { it[0].data_type == "fastq" }
        .map { meta, fastqs, _bam -> [meta, fastqs] }
        .set { ch_fastq }

    //
    // Run RNA-seq FASTQ preprocessing subworkflow
    //

    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_fastq,
        ch_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        "",
        "",
        "",
        true,
        false,
        params.skip_trimming,
        true,
        false,
        false,
        params.trimmer,
        params.min_trimmed_reads,
        params.save_trimmed,
        false,
        false,
        0,
        params.stranded_threshold,
        params.unstranded_threshold,
        false
    )

    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
    ch_versions                       = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)
    ch_trim_read_count                = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_read_count

    ch_trim_status = ch_trim_read_count
        .map {
            meta, num_reads ->
                return [ meta.id, num_reads > params.min_trimmed_reads.toFloat() ]
        }

    //
    // ALIGNMENT WITH STAR
    //

    STAR_ALIGN (
        FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads,
        ch_star_index.map { [ [:], it ] },
        ch_gtf.map { [ [:], it ] },
        false,
        "",
        ""
    )
    ch_genome_bam          = STAR_ALIGN.out.bam
    ch_transcriptome_bam   = STAR_ALIGN.out.bam_transcript
    ch_star_log            = STAR_ALIGN.out.log_final
    ch_unaligned_sequences = STAR_ALIGN.out.fastq

    ch_multiqc_files = ch_multiqc_files.mix(ch_star_log.collect{it[1]})
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    ch_star_log
        .map { meta, align_log -> [ meta ] + getStarPercentMapped(params, align_log) }
        .set { ch_percent_mapped }

    ch_map_status = ch_percent_mapped
        .map {
            meta, _mapped, pass ->
                return [ meta.id, pass ]
        }

    ch_genome_bam
        .join(ch_percent_mapped, by: [0])
        .map { meta, ofile, _mapped, pass -> if (pass) [ meta, ofile ] }
        .set { ch_genome_bam }

    ch_percent_mapped
        .branch { meta, mapped, pass ->
            pass: pass
                return [ "$meta.id\t$mapped" ]
            fail: !pass
                return [ "$meta.id\t$mapped" ]
        }
        .set { ch_pass_fail_mapped }

    ch_pass_fail_mapped
        .fail
        .collect()
        .map {
            tsv_data ->
                def header = ["Sample", "STAR uniquely mapped reads (%)"]
                sample_status_header_multiqc.text + multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_fail_mapping_multiqc }
    ch_multiqc_files = ch_multiqc_files.mix(ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv'))

    //
    // ADD BAM FILES FROM PREVIOUS ALIGNMENT DIRECTLY FROM SAMPLESHEET
    //

    ch_samplesheet
        .filter{ it[0].data_type == "bam" }
        .map {
            meta, _fastq1s, _fastq2s, bams ->
                return [ meta, bams ]
        }
        .mix ( ch_genome_bam )
        .set { ch_genome_bam }

    //
    // MODULE: Sort BAM files
    //

    SAMTOOLS_SORT ( ch_genome_bam, ch_fasta.map { [ [:], it ] } )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    //
    // SUBWORKFLOW: Mark duplicate reads
    //

    BAM_MARKDUPLICATES_PICARD (
        SAMTOOLS_SORT.out.bam,
        ch_fasta.map { [ [:], it ] },
        ch_fai.map { [ [:], it ] }
    )
    ch_genome_bam       = BAM_MARKDUPLICATES_PICARD.out.bam
    ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.bai
    ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]})

    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)

    //
    // MODULE: gene counts using featureCounts
    //

    ch_genome_bam
        .combine(ch_gtf)
        .set { ch_featurecounts }

    SUBREAD_FEATURECOUNTS_GENES(
        ch_featurecounts
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_GENES.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS_GENES.out.summary.collect{it[1]})

    GENEID_TO_GENENAME(
        ch_gtf,
        params.featurecounts_group_type,
        params.gene_name
    )
    ch_versions = ch_versions.mix(GENEID_TO_GENENAME.out.versions)

    // Check if gene names were found in GTF file, if not use NO_FILE dummy file
    gene_metadata = GENEID_TO_GENENAME.out.gene_names.ifEmpty(file("$projectDir/assets/NO_FILE", checkIfExists:true))

    SUBREAD_FEATURECOUNTS_GENES.out.counts
        .map{ it -> it[1] }
        .collect()
        .set { ch_sample_counts }

    JOIN_FEATURECOUNTS_MATRIX(
        ch_sample_counts,
        gene_metadata
    )

    //
    // MODULE: Feature biotype QC using featureCounts
    //

    ch_gtf
        .map { biotypeInGtf(it, params.featurecounts_biotype) }
        .set { biotype_in_gtf }

    // Prevent any samples from running if GTF file doesn't have a valid biotype
    ch_genome_bam
        .combine(ch_gtf)
        .combine(biotype_in_gtf)
        .filter { it[-1] }
        .map { it[0..<it.size()-1] }
        .set { ch_featurecounts }

    SUBREAD_FEATURECOUNTS_BIOTYPES (
        ch_featurecounts
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_BIOTYPES.out.versions.first())

    MULTIQC_CUSTOM_BIOTYPE (
        SUBREAD_FEATURECOUNTS_BIOTYPES.out.counts,
        ch_biotypes_header_multiqc
    )
    ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_BIOTYPE.out.tsv.collect{it[1]})
    ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())

    //
    // PICARD_COLLECTRNASEQMETRICS
    //

    PICARD_COLLECTRNASEQMETRICS(
        ch_genome_bam,
        ch_ref_flat,
        ch_fasta,
        []
    )
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTRNASEQMETRICS.out.metrics.collect{it[1]})
    ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions)

    //
    // DOWNSTREAM ANALYSIS OF COUNT MATRIX
    //

    ch_count_matrix = Channel.empty()
    if ( params.counts ) {
        ch_count_matrix = counts
        gene_column_nr = params.gene_column_nr
        gene_id_index = params.gene_id_index
    } else {
        ch_count_matrix = JOIN_FEATURECOUNTS_MATRIX.out.joined_counts
        gene_column_nr = Channel.of("7")
            .combine(GENEID_TO_GENENAME.out.gene_names)
            .map { it[0] }
            .ifEmpty("6")
            .first()
        gene_id_index = Channel.of("7")
            .combine(GENEID_TO_GENENAME.out.gene_names)
            .map { it[0] }
            .ifEmpty("1")
            .first()
    }

    if ( params.exclude_list ) {

        SAMPLE_FILTER(
            ch_count_matrix,
            params.exclude_list
        )
        ch_versions = ch_versions.mix(SAMPLE_FILTER.out.versions)

        SAMPLE_FILTER.out.filtered
            .map { file -> [["id":file.baseName], file] }
            .set { ch_counts_filt }

    } else {

        ch_count_matrix
            .map { file -> [["id":file.baseName], file] }
            .set { ch_counts_filt }

    }

    //
    // DIMENSIONALITY REDUCTION SUBWORKFLOW
    //

    DIM_REDUCTION(
        ch_counts_filt,
        gene_column_nr,
        gene_id_index,
        metadata,
        params.frac_expressed
    )
    ch_versions = ch_versions.mix(DIM_REDUCTION.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(DIM_REDUCTION.out.multiqc_files)

    // COUNT SUBSETS

    if ( params.split_file ) {

        SPLIT_COUNT_MATRIX(
            ch_counts_filt,
            gene_column_nr,
            file(params.split_file, checkIfExists: true),
        )
        ch_versions = ch_versions.mix(SPLIT_COUNT_MATRIX.out.versions)

        SPLIT_COUNT_MATRIX.out.matrices
            .flatten()
            .map { file -> [["id":file.baseName.minus(".subset")], file] }
            .set { ch_counts_split }

    } else {

        ch_counts_split = ch_counts_filt

    }

    // DEA AND FUNCTIONAL

    channel
        .fromList(params.comparisons.split(',').flatten())
        .set{ comparisons_ch }
    channel
        .fromList(params.genesets.split(',').flatten())
        .set{ genesets_ch }

    DIFFERENTIAL_EXPRESSION(
        ch_counts_split,
        gene_column_nr,
        gene_id_index,
        metadata,
        params.model_formula,
        comparisons_ch,
        genesets_ch,
        params.frac_expressed,
        params.fdr_threshold,
        params.lfc_threshold,
        params.fdr_pathways,
        params.n_pathways
    )
    ch_versions = ch_versions.mix(DIFFERENTIAL_EXPRESSION.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(DIFFERENTIAL_EXPRESSION.out.multiqc_files)

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
    ch_multiqc_report = MULTIQC.out.report

    emit:
    trim_status    = ch_trim_status    // channel: [id, boolean]
    map_status     = ch_map_status     // channel: [id, boolean]
    // strand_status  = ch_strand_status  // channel: [id, boolean]
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

//
// Function that parses and returns the alignment rate from the STAR log output
//
def getStarPercentMapped(params, align_log) {
    def percent_aligned = 0
    def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
    align_log.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_aligned = matcher[0][1].toFloat()
        }
    }

    def pass = false
    if (percent_aligned >= params.min_mapped_reads.toFloat()) {
        pass = true
    }
    return [ percent_aligned, pass ]
}

//
// Function to check whether biotype field exists in GTF file
//
def biotypeInGtf(gtf_file, biotype) {
    def hits = 0
    gtf_file.eachLine { line ->
        def attributes = line.split('\t')[-1].split()
        if (attributes.contains(biotype)) {
            hits += 1
        }
    }
    if (hits) {
        return true
    } else {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Biotype attribute '${biotype}' not found in the last column of the GTF file!\n\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  Amend '--featurecounts_group_type' to change this behaviour.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        return false
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
