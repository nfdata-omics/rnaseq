process OBJ_CONSTRUCTION {
    tag "$count_file"
    label 'process_single'

    container "docker.io/sddcunit/downstream:rnaseq-1.0.2"

    input:
    path count_file
    path metadata_table

    output:
    tuple val([id:"$count_file"]), path("*.rds"), emit: rds
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the this pipeline. Please use docker or singularity containers."
    }
    args = task.ext.args ?: ''
    """
    rnaseq_obj_setup.R $args $count_file $metadata_table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    rnaseq_obj_setup.R --versions >> versions.yml
    """

    stub:
    """
    touch rna_SummExp.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    rnaseq_obj_setup.R --versions >> versions.yml
    """
}