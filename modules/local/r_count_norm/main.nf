process R_COUNT_NORM {
    tag "$count_file"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the this pipeline. Please use docker or singularity containers."
    }
    container "docker.io/sddcunit/downstream:rnaseq-1.0.0"

    input:
    path count_file, name: 'gene_counts.tsv'

    output:
    path "lib_size_factors.cpm.txt", emit: size_factors
    path "cpm.txt"                 , emit: cpm
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    count_norm.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    """

    stub:
    """
    touch cpm.txt
    touch lib_size_factors.cpm.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    """
}
