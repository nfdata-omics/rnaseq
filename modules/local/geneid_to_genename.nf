process GENEID_TO_GENENAME {
    tag "$gtf_file.name"
    label 'process_single'

    conda "bioconda::mawk=1.3.4 conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path gtf_file
    val gtf_filter
    val gene_name

    output:
    path "GeneName_dict.txt", emit: gene_names, optional: true
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def command = gtf_file.toString().endsWith('.gz') ? 'zcat' : 'cat'
    """
#!/bin/bash

echo "Geneid\tGeneName" > GeneName_dict.txt
$command $gtf_file | grep -v "##" \
    | awk -F"\\t" '(\$3=="$gtf_filter") {print \$9}' \
    | sed -n 's/.*gene_id "\\([^ ]*\\)".*$gene_name "\\([^ ]*\\)".*/\\1\\t\\2/p' \
    | sort -u \
    >> GeneName_dict.txt

if [[ \$(wc -l < GeneName_dict.txt) -eq 1 ]]; then
    echo "No gene names found in GTF file."
    rm GeneName_dict.txt
fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    mawk: \$(mawk -W version 2> /dev/null | sed -n 's/^mawk \\([^\n]*\\).*/\\1/p')
    sed: \$(sed --version | sed -n 's/sed (GNU sed) \\([^\n]*\\).*/\\1/p')
END_VERSIONS
    """

    stub:
    """
#!/bin/bash

touch GeneName_dict.txt

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    mawk: \$(mawk -W version 2> /dev/null | sed -n 's/^mawk \\([^\n]*\\).*/\\1/p')
    sed: \$(sed --version | sed -n 's/sed (GNU sed) \\([^\n]*\\).*/\\1/p')
END_VERSIONS
    """
}
