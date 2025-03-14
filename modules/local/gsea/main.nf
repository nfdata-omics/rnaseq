process GSEA {
    tag "${meta.id}_${meta.cf}_${meta.geneset}"
    label 'process_single'

    container "docker.io/nfdata/clusterprofiler:v4.14.4"

    input:
    tuple val(meta), path(degs_table), path(gmt_file)
    val fdr_pathways

    output:
    tuple val(meta), path("gsea.*.xlsx"), emit: xlsx
    path("gsea.*.pdf"), optional: true  , emit: plot
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GSEA module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    """
    gsea.R $args --FDR $fdr_pathways --ranking logFC $degs_table $gmt_file

    echo "${task.process}": > versions.yml
    gsea.R --version >> versions.yml
    """

    stub:
    """
    touch gsea.${meta.id}_${meta.cf}.${gmt_file}.xlsx
    touch gsea.${meta.id}_${meta.cf}.${gmt_file}.pdf

    echo "${task.process}": > versions.yml
    gsea.R --version >> versions.yml
    """
}
