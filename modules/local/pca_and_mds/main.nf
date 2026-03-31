process PCA_AND_MDS {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(normalized_counts)
    val frac_samples

    output:
    path "versions.yml"                , emit: versions
    path "PCA_scores.txt"              , emit: pca
    path "PCA_explained_variance.txt"  , emit: pca_var
    path "MDS_scores.txt"              , emit: mds
    path "exprs_data_to_hc_dendrog.txt", emit: red_matrix
    path "PC1_vs_PC2_scoreplots.pdf"   , emit: PC12_pdf
    path "PC3_vs_PC4_scoreplots.pdf"   , emit: PC34_pdf
    path "MDS_scoreplots.pdf"          , emit: MDS_pdf

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PCA_AND_MDS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    """
    pca_mds.R $args -e $frac_samples --top_var 5000 $normalized_counts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    pca_mds.R --version >> versions.yml
    """

    stub:
    """
    touch PCA_scores.txt
    touch PCA_explained_variance.txt
    touch MDS_scores.txt
    touch exprs_data_to_hc_dendrog.txt
    touch PC1_vs_PC2_scoreplots.pdf
    touch PC3_vs_PC4_scoreplots.pdf
    touch MDS_scoreplots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    pca_mds.R --version >> versions.yml
    """
}
