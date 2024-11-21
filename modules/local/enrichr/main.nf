process ENRICHR {
    tag "$meta.id_$meta.cf"
    label 'process_single'

    container "docker.io/sddcunit/downstream:rnaseq-1.0.2"

    input:
    tuple val(meta), path(degs_table)
    val fdr_threshold
    val lfc_threshold

    output:
    tuple val(meta), path("*_all.xlsx") , emit: enrich_all
    tuple val(meta), path("*_up.xlsx")  , emit: enrich_up
    tuple val(meta), path("*_down.xlsx"), emit: enrich_down
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the this pipeline. Please use docker or singularity containers."
    }
    args = task.ext.args ?: ''
    """
    ora_enrichr.R $args --FDR $fdr_threshold --logfc $lfc_threshold $degs_table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    ora_enrichr.R --version >> versions.yml
    """

    stub:
    """
    touch enrichr_${meta.id}_${meta.cf}_all.xlsx
    touch enrichr_${meta.id}_${meta.cf}_up.xlsx
    touch enrichr_${meta.id}_${meta.cf}_down.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    ora_enrichr.R --version >> versions.yml
    """
}