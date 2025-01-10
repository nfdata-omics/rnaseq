
include { OBJ_CONSTRUCTION }      from '../../modules/local/obj_construction/main'
include { R_COUNT_NORM }          from '../../modules/local/r_count_norm/main'
include { CONTROL_GENE_HEATMAP }  from '../../modules/local/control_gene_heatmap/main'
include { SAMPLES_CORRELATION }   from '../../modules/local/samples_correlation/main
include { PCA_AND_MDS }           from '../../modules/local/pca_and_mds/main'

workflow DIM_REDUCTION {
    take:
    count_matrix
    gene_column_nr
    gene_id_index
    metadata_table
    frac_expressed

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

    //
    // HEATMAP OF CONTROL GENES
    //
    CONTROL_GENE_HEATMAP(
        R_COUNT_NORM.out.rds,
        "${workflow.projectDir}/assets/hsapiens_ctrl_genes.txt",
        ""
    )
    ch_versions = ch_versions.mix(CONTROL_GENE_HEATMAP.out.versions)

    //
    // CALCULATE SAMPLES CORRELATION
    //
    SAMPLES_CORRELATION(
    	R_COUNT_NORM.out.rds,
        ""
    )
    ch_versions = ch_versions.mix(SAMPLES_CORRELATION.out.versions)

    //
    // COMPUTE PCA AND MDS COORDINATES
    //
    PCA_AND_MDS(
        R_COUNT_NORM.out.rds,
        frac_expressed
    )
    ch_versions = ch_versions.mix(PCA_AND_MDS.out.versions)

    emit:
    versions    = ch_versions

}
