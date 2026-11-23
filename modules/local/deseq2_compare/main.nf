process DESEQ2_COMPARE {
    tag "${meta.id}_${meta.cf.replace('@', '.').replace('*', '_x_')}"
    label 'process_single'

    container "docker.io/nfdata/bulk_rnaseq:v1.0.1"

    input:
    tuple val(meta), path(model), val(comparison)
    val fdr_threshold

    output:
    tuple val(meta), path("deseq2_toptable.*.txt") , emit: dge
    tuple val(meta), path("deseq2_toptable.*.xlsx"), emit: dge_xlsx
    tuple val(meta), path("deseq2_summary.*.txt")  , emit: summary
    tuple val(meta), path("MAplot.*.pdf")          , emit: ma_plot
    tuple val(meta), path("pvalue_hist.*.pdf")     , emit: hist_pval
    tuple val(meta), path("volcano.*.pdf")         , emit: volcano
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DESEQ2_COMPARE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    """
    dge_deseq2_results.R $args --FDR $fdr_threshold $model "$comparison"
    for f in deseq2_toptable.*.txt; do
        comparison=\${f#deseq2_toptable.}
        comparison=\${comparison%.txt}
        cat <(printf '#DGE\\t%s\\t%s\\n' "${meta.id}" "\$comparison") "\$f" > out_tmp
        mv out_tmp "\$f"
    done

    for f in deseq2_summary.*.txt; do
        comparison=\${f#deseq2_summary.}
        comparison=\${comparison%.txt}
        cat <(printf '#DGE\\t%s\\t%s\\n' "${meta.id}" "\$comparison") "\$f" > out_tmp2
        mv out_tmp2 "\$f"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    dge_deseq2_results.R --version >> versions.yml
    """

    stub:
    """
    touch deseq2_toptable.${meta.id}_${comparison.replace('@', '.').replace('*', '_x_')}.txt
    touch deseq2_toptable.${meta.id}_${comparison.replace('@', '.').replace('*', '_x_')}.xlsx
    touch deseq2_summary.${meta.id}_${comparison.replace('@', '.').replace('*', '_x_')}.txt
    touch MAplot.${meta.id}_${comparison.replace('@', '.').replace('*', '_x_')}.pdf
    touch pvalue_hist.${meta.id}_${comparison.replace('@', '.').replace('*', '_x_')}.pdf
    touch volcano.${meta.id}_${comparison.replace('@', '.').replace('*', '_x_')}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    dge_deseq2_results.R --version >> versions.yml
    """
}
