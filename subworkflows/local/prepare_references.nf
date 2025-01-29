
include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_REFFLAT          } from '../../modules/nf-core/gunzip'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../modules/nf-core/untar'

include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from '../../modules/nf-core/rsem/preparereference'
include { CUSTOM_GETCHROMSIZES              } from '../../modules/nf-core/custom/getchromsizes'
include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/star/genomegenerate'
include { SALMON_INDEX                      } from '../../modules/nf-core/salmon/index'
include { UCSC_GTFTOGENEPRED                } from '../../modules/nf-core/ucsc/gtftogenepred/main'

include { GTF_FILTER                        } from '../../modules/local/gtf_filter'
include { GTF2BED                           } from '../../modules/local/gtf2bed'
workflow PREPARE_REFERENCES {
    take:
    fasta                    //      file: /path/to/genome.fasta
    gtf                      //      file: /path/to/genome.gtf
    transcript_fasta         //      file: /path/to/transcript.fasta
    star_index               // directory: /path/to/star/index/
    salmon_index             // directory: /path/to/salmon/index/
    refflat                  //      file: /path/to/refFlat.txt

    main:
    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], file(fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta, checkIfExists: true))
    }

    //
    // Uncompress GTF annotation file if required
    //
    if (gtf.endsWith('.gz')) {
        ch_gtf      = GUNZIP_GTF ( [ [:], file(gtf, checkIfExists: true) ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = Channel.value(file(gtf, checkIfExists: true))
    }

    //
    // Filter GTF annotation file
    //
    GTF_FILTER ( ch_fasta, ch_gtf )
    ch_gtf = GTF_FILTER.out.genome_gtf
    ch_versions = ch_versions.mix(GTF_FILTER.out.versions)

    //
    // Create gene BED annotation file
    //
    ch_gene_bed = GTF2BED ( ch_gtf ).bed
    ch_versions = ch_versions.mix(GTF2BED.out.versions)

    //
    // Uncompress transcript fasta file / create if required
    //
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], file(transcript_fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta, checkIfExists: true))
        }
    } else {
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Uncompress STAR index or generate from scratch if required
    //
    if (star_index) {
        if (star_index.endsWith('.tar.gz')) {
            ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map { it[1] }
            ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
        } else {
            ch_star_index = Channel.value(file(star_index))
        }
    } else {
        ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] } ).index.map { it[1] }
        ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    if (salmon_index) {
        if (salmon_index.endsWith('.tar.gz')) {
            ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], salmon_index ] ).untar.map { it[1] }
            ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
        } else {
            ch_salmon_index = Channel.value(file(salmon_index))
        }
    } else {
        ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
        ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
    }

    //
    // Uncompress refFlat file or generate it from the gtf
    //
    if (refflat) {
        if (refflat.endsWith('.gz')) {
            ch_refflat = GUNZIP_REFFLAT ( [ [:], file(refflat, checkIfExists: true) ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_REFFLAT.out.versions)
        } else {
            ch_refflat = Channel.value(file(refflat, checkIfExists: true))
        }
    } else {
        ch_refflat = UCSC_GTFTOGENEPRED ( ch_gtf.map { [ [:], it ] } ).refflat.map { it[1] }
        ch_versions = ch_versions.mix(UCSC_GTFTOGENEPRED.out.versions)
    }

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    fai              = ch_fai                    // channel: path(genome.fai)
    gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    star_index       = ch_star_index             // channel: path(star/index/)
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)
    refflat          = ch_refflat                // channel: path(genome.refFlat)
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

}