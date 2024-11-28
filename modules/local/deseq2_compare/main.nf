process DESEQ2_COMPARE {
    tag "${meta.id}_$comparison"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(model)
    val comparison

    output:
    tuple val(meta), path("deseq2_toptable.*.txt"), emit: dge
    path "deseq2_summary.*.txt"                   , emit: summary
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the this pipeline. Please use docker or singularity containers."
    }
    args = task.ext.args ?: ''
    """
    dge_deseq2_results.R $args $model $comparison

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    dge_deseq2_results.R --version >> versions.yml
    """

    stub:
    """
    touch deseq2_toptable.${meta.id}_${comparison}.txt
    touch deseq2_summary.${meta.id}_${comparison}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    dge_deseq2_results.R --version >> versions.yml
    """
}