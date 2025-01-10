process SAMPLES_CORRELATION {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(normalized_counts)
    val covariates_highlight

    output:
    path "samples_correlation_heatmap.pdf", emit: corr_heatmap_pdf
    path "samples_correlation_table.txt"  , emit: corr_matrix
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the this pipeline. Please use docker or singularity containers."
    }
    args = task.ext.args ?: ''
    cov_list = covariates_highlight ? "-a $covariates_highlight" : ''
    """
    samples_correlation.R $args $cov_list $normalized_counts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    samples_correlation.R --version >> versions.yml
    """

    stub:
    """
    touch samples_correlation_table.txt
    touch samples_correlation_heatmap.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    samples_correlation.R --version >> versions.yml
    """
}
