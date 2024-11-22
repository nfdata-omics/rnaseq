process R_COUNT_NORM {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/sddcunit/downstream:rnaseq-1.0.2"

    input:
    tuple val(meta), path(annotated_counts)

    output:
    tuple val(meta), path("*.norm.rds")     , emit: rds
    path "lib_size_factors.txt"             , emit: size_factors
    path "cpm.txt"                          , emit: cpm
    path 'log.norm.rpkm.txt', optional: true, emit: rpkm
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the this pipeline. Please use docker or singularity containers."
    }
    args = task.ext.args ?: ''
    """
    count_norm.R $args $annotated_counts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    count_norm.R --version >> versions.yml
    """

    stub:
    """
    touch cpm.txt
    touch lib_size_factors.cpm.txt
    touch rnaseq_matrix.norm.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    count_norm.R --version >> versions.yml
    """
}
