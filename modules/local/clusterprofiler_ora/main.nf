process CLUSTERPROFILER_ORA {
    tag "${meta.id}_${meta.cf}_${meta.geneset}"
    label 'process_single'

    container "docker.io/nfdata/clusterprofiler:v4.14.4"

    input:
    tuple val(meta), path(degs_table), path(gmt_file)
    val fdr_threshold
    val lfc_threshold

    output:
    tuple val(meta), path("ora_CP*.all.txt") , optional: false, emit: CP_all
    tuple val(meta), path("ora_CP*.up.txt")  , optional: false, emit: CP_up
    tuple val(meta), path("ora_CP*.down.txt"), optional: false, emit: CP_down
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CLUSTERPROFILER_ORA module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    """
    ora_clusterProfiler.R $args --FDR $fdr_threshold --logfc $lfc_threshold $degs_table $gmt_file

    echo "${task.process}": > versions.yml
    ora_clusterProfiler.R --version >> versions.yml
    """

    stub:
    """
    touch ora_CP.${meta.id}_${meta.cf}.${gmt_file}.all.txt
    touch ora_CP.${meta.id}_${meta.cf}.${gmt_file}.up.txt
    touch ora_CP.${meta.id}_${meta.cf}.${gmt_file}.down.txt

    echo "${task.process}": > versions.yml
    ora_clusterProfiler.R --version >> versions.yml
    """
}
