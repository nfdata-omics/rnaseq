process JOIN_FEATURECOUNTS_MATRIX {
    tag "joined_matrix"
    label 'process_single'

    conda "conda-forge::pandas==2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path count_matrix
    path gene_metadata

    output:
    path "joined_matrix.tsv", emit: joined_counts
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genename_arg = gene_metadata.name != 'NO_FILE' ? "--genemeta $gene_metadata" : ''
    """
#!/bin/bash

join_featurecounts_matrix.py \
    ${genename_arg} --output joined_matrix.tsv \
    ${count_matrix}

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | awk '{print \$2}')
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
END_VERSIONS
    """

    stub:
    """
#!/bin/bash

touch joined_matrix.tsv

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | awk '{print \$2}')
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
END_VERSIONS
    """
}
