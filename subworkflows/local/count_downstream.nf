
include { OBJ_CONSTRUCTION }      from '../../modules/local/obj_construction/main'
include { R_COUNT_NORM }          from '../../modules/local/r_count_norm/main'
include { CONTROL_GENE_HEATMAP }  from '../../modules/local/control_gene_heatmap/main'
include { DIM_REDUCTION }         from '../../modules/local/dim_reduction/main'
include { DESEQ2_FIT }            from '../../modules/local/deseq2_fit/main'
include { DESEQ2_COMPARE }        from '../../modules/local/deseq2_compare/main'
include { ENRICHR }               from '../../modules/local/enrichr/main'

workflow COUNT_DOWNSTREAM {
    take:
    count_matrix
    metadata_table
    model_formula
    comparisons_ch
    frac_expressed
    fdr_threshold
    lfc_threshold

    main:

    //
    // R-OBJECT CONSTRUCTION
    //
    // NOTE: metadata needs to be manually curated before execution: numerical variables are continous covariates
    // while string var are categorical
    //
    OBJ_CONSTRUCTION(
        count_matrix,
        metadata_table
    )

    //
    // COUNT NORMALIZATION
    //
    R_COUNT_NORM(
        OBJ_CONSTRUCTION.out.rds
    )

    //
    // HEATMAP OF CONTROL GENES
    //
    CONTROL_GENE_HEATMAP(
        R_COUNT_NORM.out.rds,
        "${workflow.projectDir}/assets/hsapiens_ctrl_genes.txt",
        ""
    )

    //
    // DIMENSIONALITY REDUCTION
    //
    DIM_REDUCTION(
        R_COUNT_NORM.out.rds,
        frac_expressed
    )

    //
    // DESEQ2 MODEL FIT
    //
    DESEQ2_FIT(
        R_COUNT_NORM.out.rds,
        model_formula,
        frac_expressed
    )

    //
    // DESEQ2 COMPARISONS
    //
    DESEQ2_COMPARE(
        DESEQ2_FIT.out.rds,
        comparisons_ch
    )

    //
    // FUNCTIONAL ANALYSIS: ENRICHR
    //
    ENRICHR(
        DESEQ2_COMPARE.out.dge,
        fdr_threshold,
        lfc_threshold
    )

}