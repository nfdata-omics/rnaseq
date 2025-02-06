
include { OBJ_CONSTRUCTION }      from '../../modules/local/obj_construction/main'
include { R_COUNT_NORM }          from '../../modules/local/r_count_norm/main'
include { DESEQ2_FIT }            from '../../modules/local/deseq2_fit/main'
include { DESEQ2_COMPARE }        from '../../modules/local/deseq2_compare/main'
include { ENRICHR }               from '../../modules/local/enrichr/main'
include { ENRICHR_TOPN }          from '../../modules/local/enrichr_topn/main'

workflow DIFFERENTIAL_EXPRESSION {
    take:
    count_matrix
    gene_column_nr
    gene_id_index
    metadata_table
    model_formula
    comparisons_ch
    frac_expressed
    fdr_threshold
    lfc_threshold
    fdr_pathways
    n_pathways

    main:

    ch_versions = Channel.empty()

    //
    // R-OBJECT CONSTRUCTION
    //
    // NOTE: metadata needs to be manually curated before execution: numerical variables are continous covariates
    // while string var are categorical
    //
    OBJ_CONSTRUCTION(
        count_matrix,
        gene_column_nr,
        gene_id_index,
        metadata_table
    )
    ch_versions = ch_versions.mix(OBJ_CONSTRUCTION.out.versions)

    //
    // COUNT NORMALIZATION
    //
    R_COUNT_NORM(
        OBJ_CONSTRUCTION.out.rds
    )
    ch_versions = ch_versions.mix(R_COUNT_NORM.out.versions)

    if ( model_formula ) {

        //
        // DESEQ2 MODEL FIT
        //
        DESEQ2_FIT(
            R_COUNT_NORM.out.rds,
            model_formula,
            frac_expressed
        )
        ch_versions = ch_versions.mix(DESEQ2_FIT.out.versions)

        if ( comparisons_ch ) {

            DESEQ2_FIT.out.rds
                .combine(comparisons_ch)
                .map{ meta, rds, cf -> [["id": meta.id, "cf": cf], rds, cf] }
                .set{ ch_rds_cf }

            //
            // DESEQ2 COMPARISONS
            //
            DESEQ2_COMPARE(
                ch_rds_cf,
                fdr_threshold
            )
            ch_versions = ch_versions.mix(DESEQ2_COMPARE.out.versions)

            //
            // FUNCTIONAL ANALYSIS: ENRICHR
            //
            ENRICHR(
                DESEQ2_COMPARE.out.dge,
                fdr_threshold,
                lfc_threshold
            )
            ch_versions = ch_versions.mix(ENRICHR.out.versions)

            ENRICHR.out.enrich_all
                .mix(
                    ENRICHR.out.enrich_up,
                    ENRICHR.out.enrich_down
                )
                .set{ ch_enrichr }

            //
            // ENRICHR TOP RESULTS EXTRACTION
            //
            ENRICHR_TOPN(
                ch_enrichr,
                fdr_pathways,
                n_pathways
            )
            ch_versions = ch_versions.mix(ENRICHR_TOPN.out.versions)

        }
    }

    emit:
    versions    = ch_versions
}
