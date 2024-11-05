process CONTROL_GENE_HEATMAP {
    tag "$count_file"
    label 'process_single'

    container "docker.io/sddcunit/downstream:rnaseq-1.0.1"

    input:
    path count_file
    path gene_list
    path metadata_table
    val covariates_highlight

    output:
    path "heatmap_control_genes.pdf",  emit: heatmap_pdf
    path "countrol_genes_exprs.txt",   emit: count_subset
    path "versions.yml"            ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the this pipeline. Please use docker or singularity containers."
    }
    args = task.ext.args ?: ''
    """
    control_genes.R $args -a $covariates_highlight \
        $count_file $gene_list $metadata_table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    control_genes.R --versions >> versions.yml
    """

    stub:
    """
    touch countrol_genes_exprs.txt
    touch heatmap_control_genes.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    control_genes.R --versions >> versions.yml
    """
}