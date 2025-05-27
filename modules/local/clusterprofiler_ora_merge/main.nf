process CLUSTERPROFILER_ORA_MERGE {
    tag "${meta.id}_${meta.cf}"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(list_of_tables)
    val fdr_pathways
    val n_pathways

    output:
    tuple val(meta), path("ora_CP.*.xlsx") , optional: true, emit: CP_merged
    tuple val(meta), path("ora_CP.*.txt")  , optional: true, emit: CP_topn
    tuple val(meta), path("ora_CP.*.pdf")  , optional: true, emit: CP_topn_dotplot
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CLUSTERPROFILER_ORA_MERGE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    """
    ora_clusterProfiler_merge.R $args --FDR $fdr_pathways --num $n_pathways $list_of_tables

    echo "${task.process}": > versions.yml
    ora_clusterProfiler_merge.R --version >> versions.yml
    """

    stub:
    """
    touch ora_CP.${meta.id}_${meta.cf}.TOP_N.txt
    touch ora_CP.${meta.id}_${meta.cf}.xlsx
    touch ora_CP.${meta.id}_${meta.cf}.TOP_N.pdf

    echo "${task.process}": > versions.yml
    ora_clusterProfiler_merge.R --version >> versions.yml
    """
}


