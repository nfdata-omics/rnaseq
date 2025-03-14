process ENRICHR_TOPN {
    tag "${meta.id}_${meta.cf}"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(enrich_table)
    val fdr_pathways
    val n_pathways

    output:
    tuple val(meta), path("*_TOP*.txt") , optional: true, emit: enrich_topn
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "ENRICHR_TOPN module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    """
    ora_top_results.R $args --FDR $fdr_pathways --num $n_pathways $enrich_table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    ora_top_results.R --version >> versions.yml
    """

    stub:
    """
    touch ${enrich_table.baseName}_TOP${n_pathways}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    ora_top_results.R --version >> versions.yml
    """
}
