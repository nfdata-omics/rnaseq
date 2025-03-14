process SAMPLES_CORRELATION {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(normalized_counts)
    val covariates_highlight

    output:
    path "samples_correlation_heatmap.pdf"    , emit: corr_heatmap_pdf
    path "samples_correlation_heatmap_mqc.png", emit: corr_heatmap_png
    path "samples_correlation_table.txt"      , emit: corr_matrix
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SAMPLES_CORRELATION module does not support Conda. Please use Docker / Singularity / Podman instead."
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
    touch samples_correlation_heatmap_mqc.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    samples_correlation.R --version >> versions.yml
    """
}
