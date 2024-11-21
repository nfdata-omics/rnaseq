process SPLIT_COUNT_MATRIX {
    tag "$count_file"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::pandas==2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path count_file
    val gene_column_nr
    path metadata
    val split_variable

    output:
    path "*.subset.tsv", emit: matrices
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import sys
    import os
    import pandas as pd

    count_file = "$count_file"
    gene_column_nr = $gene_column_nr
    metadata = "$metadata"
    split_variable = "$split_variable"

    # read count matrix
    df = pd.read_csv(count_file, sep='\t')
    samples = df.columns[gene_column_nr:].to_list()
    print(f"Read count matrix with {df.shape[0]} rows and {df.shape[1]} columns.")
    print(f"The matrix contains {len(samples)} sample columns.\n")

    # identify samples present in the count matrix
    gene_metadata = df.columns[:gene_column_nr].to_list()
    print(f"Count annotation columns: {gene_metadata}\nThese columns will be maintained in all the subset matrix files\n")

    # read and filter metadata table
    df_meta = pd.read_csv(metadata)
    print(f"Read metadata table with {df_meta.shape[0]} rows and {df_meta.shape[1]} columns.")
    df_meta = df_meta[df_meta.iloc[:, 0].isin(samples)]
    print(f"{len(df_meta)} samples found in the metadata table.\n")

    # extract matrix subset corresponding to each value of the split_variable metadata and write it to disk
    for meta_val in df_meta[split_variable].unique():
        print(f"Splitting full count matrix for {split_variable} = {meta_val}")
        sample_subset = df_meta[df_meta[split_variable] == meta_val].iloc[:, 0].to_list()
        count_matrix_subset = df[gene_metadata + sample_subset]
        print(f"Subset matrix has {count_matrix_subset.shape[0]} rows and {count_matrix_subset.shape[1]} columns.")
        out_fname = f"{split_variable}_{meta_val}.subset.tsv"
        count_matrix_subset.to_csv(out_fname, sep="\t", header=True, index=False)
        print(f"Matrix written to file {out_fname}.\n")

    # write version info to yaml file
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":')
        f.write(f"   python: {".".join([str(x) for x in sys.version_info[0:3]])}")
        f.write(f"   pandas: {pd.__version__}")

    """

    stub:
    """
    #!/usr/bin/env python3

    import sys
    import os
    import pandas as pd

    # Creates an empty count matrix file
    out_name = os.path.splitext("$count_file")[0] + ".subset.tsv"
    with open(out_name, 'w') as fp:
        pass

    # write version info
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":')
        f.write(f"   python: {".".join([str(x) for x in sys.version_info[0:3]])}")
        f.write(f"   pandas: {pd.__version__}")

    """
}