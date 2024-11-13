process DESEQ2_FIT {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/sddcunit/downstream:rnaseq-1.0.2"

    input:
    tuple val(meta), path(counts)
    val model_formula
    val frac_samples

    output:
    tuple val(meta), path("deseq2_obj.rds"), emit: rds
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the this pipeline. Please use docker or singularity containers."
    }
    args = task.ext.args ?: ''
    """
    dge_deseq2_fit.R $args -e $frac_samples $counts $model_formula

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    dge_deseq2_fit.R --versions >> versions.yml
    """

    stub:
    """
    touch deseq2_obj.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    dge_deseq2_fit.R --versions >> versions.yml
    """
}