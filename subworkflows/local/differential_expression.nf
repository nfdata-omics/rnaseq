
include { OBJ_CONSTRUCTION }          from '../../modules/local/obj_construction/main'
include { R_COUNT_NORM }              from '../../modules/local/r_count_norm/main'
include { DESEQ2_FIT }                from '../../modules/local/deseq2_fit/main'
include { DESEQ2_COMPARE }            from '../../modules/local/deseq2_compare/main'
//include { ENRICHR }                   from '../../modules/local/enrichr/main'
//include { ENRICHR_TOPN }              from '../../modules/local/enrichr_topn/main'
include { GSEA }                      from '../../modules/local/gsea/main'
include { GSEA_MERGE }                from '../../modules/local/gsea_merge/main'
include { CLUSTERPROFILER_ORA }       from '../../modules/local/clusterprofiler_ora/main'
include { CLUSTERPROFILER_ORA_MERGE } from '../../modules/local/clusterprofiler_ora_merge/main'
include { SUMMARY_TABLE }             from '../../modules/local/summary_table'

workflow DIFFERENTIAL_EXPRESSION {
    take:
    count_matrix
    gene_column_nr
    gene_id_index
    metadata_table
    model_formula
    ch_comparisons
    ch_genesets
    frac_expressed
    fdr_threshold
    lfc_threshold
    fdr_pathways
    n_pathways

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

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

        DESEQ2_FIT.out.rds
            .combine(ch_comparisons)
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
        ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_COMPARE.out.dge.collect{ _meta, file -> file })
        ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_COMPARE.out.summary.collect{ _meta, file -> file })

        DESEQ2_COMPARE.out.summary
            .map{ meta, result -> [meta.id, meta.cf, result] }
            .groupTuple()
            .map{ key, cfs, results -> [["id": key], cfs, results] }
            .set{ ch_dge_summary }

        SUMMARY_TABLE(
            ch_dge_summary,
        )

        DESEQ2_COMPARE.out.dge
        .combine(ch_genesets)
        .map{ meta, dge, geneset -> [ meta + ["geneset": file(geneset).baseName], dge, geneset] }
        .set{ ch_dge_geneset }

        //
        // FUNCTIONAL ANALYSIS: GSEA
        //
        GSEA(
            ch_dge_geneset,
            fdr_pathways
        )
        ch_versions = ch_versions.mix(GSEA.out.versions)

        GSEA.out.xlsx
        .map{ meta, result -> [[meta.id, meta.cf], meta, result] }
        .groupTuple()
        .map{ _key, metas, results -> [["id": metas[0].id, "cf": metas[0].cf], results] }
        .set{ ch_gsea }

        //
        // GSEA MERGING
        //
        GSEA_MERGE(
            ch_gsea
        )
        ch_versions = ch_versions.mix(GSEA_MERGE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(GSEA_MERGE.out.gsea_topn.collect{ _meta, file -> file })

        //
        // FUNCTIONAL ANALYSIS: CLUSTERPROFILER OVER-REPRESENTATION
        //
        CLUSTERPROFILER_ORA(
            ch_dge_geneset,
            fdr_threshold,
            lfc_threshold
        )
        ch_versions = ch_versions.mix(CLUSTERPROFILER_ORA.out.versions)

        CLUSTERPROFILER_ORA.out.CP_all
            .map{ meta, result -> [[meta.id, meta.cf], meta, result] }
            .groupTuple()
            .map{ _key, metas, results -> [["id": metas[0].id, "cf": metas[0].cf], results] }
            .set{ ch_CP_all }
        CLUSTERPROFILER_ORA.out.CP_up
            .map{ meta, result -> [[meta.id, meta.cf], meta, result] }
            .groupTuple()
            .map{ _key, metas, results -> [["id": metas[0].id, "cf": metas[0].cf], results] }
            .set{ ch_CP_up }
        CLUSTERPROFILER_ORA.out.CP_down
            .map{ meta, result -> [[meta.id, meta.cf], meta, result] }
            .groupTuple()
            .map{ _key, metas, results -> [["id": metas[0].id, "cf": metas[0].cf], results] }
            .set{ ch_CP_down }

        ch_CP_all
            .mix(
                ch_CP_up,
                ch_CP_down
            )
            .set{ ch_CP }

        //
        // CLUSTERPROFILER OVER-REPRESENTATION MERGING
        //
        CLUSTERPROFILER_ORA_MERGE(
            ch_CP,
            fdr_pathways,
            n_pathways
        )
        ch_versions = ch_versions.mix(CLUSTERPROFILER_ORA.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(CLUSTERPROFILER_ORA_MERGE.out.CP_topn.collect{it[1]})
    }

            //
            // FUNCTIONAL ANALYSIS: ENRICHR
            //
            /*
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
            */

    emit:
    versions    = ch_versions
    multiqc_files   = ch_multiqc_files
}
