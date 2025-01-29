process DESEQ2_FIT {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

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
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DESEQ2_FIT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    """
    dge_deseq2_fit.R $args -e $frac_samples $counts \"$model_formula\"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    dge_deseq2_fit.R --version >> versions.yml
    """

    stub:
    """
    touch deseq2_obj.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    dge_deseq2_fit.R --version >> versions.yml
    """
}