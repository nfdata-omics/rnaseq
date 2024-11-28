process CONTROL_GENE_HEATMAP {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(normalized_counts)
    path gene_list
    val covariates_highlight

    output:
    path "heatmap_control_genes.pdf", emit: heatmap_pdf
    path "control_genes_exprs.txt",   emit: count_subset
    path "versions.yml"            ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the this pipeline. Please use docker or singularity containers."
    }
    args = task.ext.args ?: ''
    cov_list = covariates_highlight ? "-a $covariates_highlight" : ''
    """
    control_genes.R $args -s $cov_list \
        $normalized_counts $gene_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    control_genes.R --version >> versions.yml
    """

    stub:
    """
    touch control_genes_exprs.txt
    touch heatmap_control_genes.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    control_genes.R --version >> versions.yml
    """
}