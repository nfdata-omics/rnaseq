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
import csv
import pandas as pd

count_file = "$count_file"
gene_column_nr = $gene_column_nr
split_file = "$split_file"

# read count matrix
df = pd.read_csv(count_file, sep='\\t')

# identify samples present in the count matrix
samples = df.columns[gene_column_nr:].to_list()
print(f"Read count matrix with {df.shape[0]} rows and {df.shape[1]} columns.")
print(f"The matrix contains {len(samples)} sample columns.\\n")

# identify gene information present in the count matrix
gene_metadata = df.columns[:gene_column_nr].to_list()
print(f"Count annotation columns: {gene_metadata}\\nThese columns will be maintained in all the subset matrix files\\n")

with open(split_file, "r", newline="") as fh:
    reader = csv.reader(fh)
    for line_nr, row in enumerate(reader, start=1):
        # Skip empty/blank lines
        if not row or all((x.strip() == "") for x in row):
            continue
        # Name of the sample subset in the first field
        subset_name = row[0].strip()
        if subset_name == "":
            print(f"Line {line_nr}: empty subset name, skipping.")
            continue
        # Remaining fields are samples
        sample_subset = [x.strip() for x in row[1:] if x.strip() != ""]
        print(f"Subset {subset_name} of the full count matrix")

        # Evaluate overlap
        actual_samples = set(sample_subset) & set(samples)
        missing_samples = set(sample_subset) - set(samples)
        print(f"  - Found {len(actual_samples)} matching samples in the count matrix.")
        if missing_samples:
            print(f"  - Warning: {len(missing_samples)} samples in the list are not present in the count matrix.")
            print(f"    Missing samples: {', '.join(missing_samples)}")

        # Subset count matrix with the available samples
        count_matrix_subset = df[gene_metadata + list(actual_samples)]
        print(f"Subset matrix has {count_matrix_subset.shape[0]} rows and {count_matrix_subset.shape[1]} columns.")
        out_fname = f"{subset_name}.subset.tsv"
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
import csv
import pandas as pd

split_file = "$split_file"

with open(split_file, "r", newline="") as fh:
    reader = csv.reader(fh)
    for line_nr, row in enumerate(reader, start=1):
        # Skip empty/blank lines
        if not row or all((x.strip() == "") for x in row):
            continue
        # Name of the sample subset in the first field
        subset_name = row[0].strip()
        if subset_name == "":
            print(f"Line {line_nr}: empty subset name, skipping.")
            continue
        # Creating empty subset matrix
        out_fname = f"{subset_name}.subset.tsv"
        with open(out_fname, 'w') as fp:
            pass

# write version info
with open("versions.yml", "w") as f:
    f.write('"${task.process}":')
    f.write(f"   python: {".".join([str(x) for x in sys.version_info[0:3]])}")
    f.write(f"   pandas: {pd.__version__}")


    """
}
