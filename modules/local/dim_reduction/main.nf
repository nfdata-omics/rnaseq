process DIM_REDUCTION {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/sddcunit/downstream:rnaseq-1.0.2"

    input:
    tuple val(meta), path(normalized_counts)
    val frac_samples

    output:
    path "versions.yml"                , emit: versions
    path "PCA_scores.txt"              , emit: pca
    path "PCA_explained_variance.txt"  , emit: pca_var
    path "MDS_scores.txt"              , emit: mds
    path "exprs_data_to_hc_dendrog.txt", emit: red_matrix

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the this pipeline. Please use docker or singularity containers."
    }
    args = task.ext.args ?: ''
    """
    pca_mds.R $args -e $frac_samples --top_var 5000 $normalized_counts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    pca_mds.R --versions >> versions.yml
    """

    stub:
    """
    touch PCA_scores.txt
    touch PCA_explained_variance.txt
    touch MDS_scores.txt
    touch exprs_data_to_hc_dendrog.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    pca_mds.R --versions >> versions.yml
    """
}