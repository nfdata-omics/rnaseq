process SAMPLE_FILTER {
    tag "$count_file"
    label 'process_single'

    conda "conda-forge::pandas==2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path count_file
    path samples_to_exclude

    output:
    path "*.filtered.tsv", emit: filtered
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
#!/usr/bin/env python3

import os
import sys
import pandas as pd

# read count matrix
df = pd.read_csv("$count_file", sep="\\t")
print(f"Read count matrix with {df.shape[0]} rows and {df.shape[1]} columns.\\n")

# read from input file the list of samples that should be excluded from the count matrix
with open("$samples_to_exclude") as f:
    s_exclude = list(filter(None, f.read().split("\\n")))
print(f"Filtering out samples from this list:\\n{"\\n".join(s_exclude)}\\n")

# filter out columns from s_exclude
df = df.loc[:, (~ df.columns.isin(s_exclude))]
print(f"Output count matrix with {df.shape[0]} rows and {df.shape[1]} columns.\\n")

# write matrix to output
out_name = os.path.splitext("$count_file")[0] + ".filtered.tsv"
df.to_csv(out_name, sep='\\t', header=True, index=False)

# write version info
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f"   python: {".".join([str(x) for x in sys.version_info[0:3]])}\\n")
    f.write(f"   pandas: {pd.__version__}\\n")

    """

    stub:
    """
#!/usr/bin/env python3

# Creates an empty count matrix file
out_name = os.path.splitext("$count_file")[0] + ".filtered.tsv"
with open(out_name, 'w') as fp:
    pass

# write version info
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f"   python: {".".join([str(x) for x in sys.version_info[0:3]])}\\n")
    f.write(f"   pandas: {pd.__version__}\\n")

    """
}
