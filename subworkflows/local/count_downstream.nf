
include { R_COUNT_NORM }          from '../../modules/local/r_count_norm/main'
include { CONTROL_GENE_HEATMAP }  from '../../modules/local/control_gene_heatmap/main'

workflow COUNT_DOWNSTREAM {
    take:
    count_matrix

    main:

    //
    // COUNT NORMALIZATION
    //
    R_COUNT_NORM(
        count_matrix
    )

    //
    // HEATMAP OF CONTROL GENES
    //
    CONTROL_GENE_HEATMAP(
        R_COUNT_NORM.out.cpm
    )


}