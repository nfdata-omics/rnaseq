process SPLIT_COUNT_MATRIX {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas==2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(count_file)
    val gene_column_nr
    path split_file

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
split_file = "$split_file"

# read count matrix
df = pd.read_csv(count_file, sep='\\t')
samples = df.columns[gene_column_nr:].to_list()
print(f"Read count matrix with {df.shape[0]} rows and {df.shape[1]} columns.")
print(f"The matrix contains {len(samples)} sample columns.\\n")

# identify samples present in the count matrix
gene_metadata = df.columns[:gene_column_nr].to_list()
print(f"Count annotation columns: {gene_metadata}\\nThese columns will be maintained in all the subset matrix files\\n")

# read file with lists of samples to extract to subset the count matrix
sample_lists = pd.read_csv(split_file, index_col=0, header=None).T

# loop over sample lists and subset matrix
for c in sample_lists.columns:
    print(f"Subset {c} of the full count matrix")
    sample_subset = sample_lists[c].dropna().tolist()
    count_matrix_subset = df[gene_metadata + sample_subset]
    print(f"Subset matrix has {count_matrix_subset.shape[0]} rows and {count_matrix_subset.shape[1]} columns.")
    out_fname = f"{c}.subset.tsv"
    count_matrix_subset.to_csv(out_fname, sep="\\t", header=True, index=False)
    print(f"Matrix written to file {out_fname}.\\n")

# write version info to yaml file
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f"   python: {".".join([str(x) for x in sys.version_info[0:3]])}\\n")
    f.write(f"   pandas: {pd.__version__}\\n")

    """

    stub:
    """
    #!/usr/bin/env python3

    import sys
    import os
    import pandas as pd

    split_file = "$split_file"

    # read file with lists of samples to extract to subset the count matrix
    sample_lists = pd.read_csv(split_file, index_col=0, header=None).T

    # loop over sample lists and create empty subset matrix
    for c in sample_lists.columns:
        out_fname = f"{c}.subset.tsv"
        with open(out_fname, 'w') as fp:
            pass

    # write version info
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":')
        f.write(f"   python: {".".join([str(x) for x in sys.version_info[0:3]])}")
        f.write(f"   pandas: {pd.__version__}")

    """
}