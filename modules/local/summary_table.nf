process SUMMARY_TABLE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), val(cf_names), path(single_cf_summaries)

    output:
    path "dea_*_summary.txt", emit: joint_summary
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Example: Loop over two lists (values and files) and print them side by side
    names=(${cf_names.join(" ")})
    files=(${single_cf_summaries})

    echo -e 'Comparison\tUpreg\tUpreg_percent\tDownreg\tDownreg_percent' > "dea_${meta.id}_summary.txt"

    for idx in "\${!names[@]}"; do
        comparison=\$(echo "\${names[\$idx]}" | sed -E 's/\\//_/g; s/^([^:]+):([^:]+):([^:]+)\$/\\1_\\2_vs_\\3/')
        upreg=\$(cat \${files[\$idx]} | grep up | sed 's/.*: //' | tr -d ' ' | tr ',' '\t')
        downreg=\$(cat \${files[\$idx]} | grep down | sed 's/.*: //' | tr -d ' ' | tr ',' '\t')
        echo -e "\${comparison}\t\${upreg}\t\${downreg}" >> "dea_${meta.id}_summary.txt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version | sed -n 's/sed (GNU sed) \\([^\\n]*\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    """
    touch "dea_${meta.id}_summary.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version | sed -n 's/sed (GNU sed) \\([^\\n]*\\).*/\\1/p')
    END_VERSIONS
    """
}
