process OBJ_CONSTRUCTION {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(count_file)
    val gene_column_nr
    val gene_id_index
    path metadata_table

    output:
    tuple val(meta), path("*.rds"), emit: rds
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "OBJ_CONSTRUCTION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    """
    rnaseq_obj_setup.R $args -i $gene_id_index -c $gene_column_nr $count_file $metadata_table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    rnaseq_obj_setup.R --version >> versions.yml
    """

    stub:
    """
    touch rna_SummExp.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    rnaseq_obj_setup.R --version >> versions.yml
    """
}