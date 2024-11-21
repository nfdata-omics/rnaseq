process DESEQ2_COMPARE {
    tag "$meta.id_$comparison"
    label 'process_single'

    container "docker.io/sddcunit/downstream:rnaseq-1.0.2"

    input:
    tuple val(meta), path(model)
    val comparison

    output:
    tuple val(meta), path("*.txt"), emit: dge
    path "versions.yml"           , emit: versions

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
    touch dge.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    dge_deseq2_results.R --version >> versions.yml
    """
}