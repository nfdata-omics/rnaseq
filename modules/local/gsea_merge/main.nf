process GSEA_MERGE {
    tag "${meta.id}_${meta.cf}"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(list_of_tables)

    output:
    tuple val(meta), path("gsea*.xlsx")  , emit: gsea_merged
    tuple val(meta), path("gsea*.txt")   , emit: gsea_topn
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GSEA_MERGE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    """
    gsea_merge.R $args $list_of_tables

    echo "${task.process}": > versions.yml
    gsea_merge.R --version >> versions.yml
    """

    stub:
    """
    touch gsea.${meta.id}_${meta.cf}.xlsx
    touch gsea.${meta.id}_${meta.cf}.TOP_N.txt

    echo "${task.process}": > versions.yml
    gsea_merge.R --version >> versions.yml
    """
}
