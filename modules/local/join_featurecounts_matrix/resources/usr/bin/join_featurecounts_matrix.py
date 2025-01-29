#!/usr/bin/env python

import argparse
import pandas as pd

def join_featurecounts_matrix(count_files, output_file, genemeta_file=None):
    """
    Joins multiple featureCounts output files into a single matrix and optionally merges
    it with gene metadata.

    Parameters:
    count_files (list of str): List of file paths to the featureCounts output files to be joined.
    output_file (str): File path where the merged output will be saved.
    genemeta_file (str, optional): File path to the gene metadata file to be merged with the count
    data. Default is None.

    Returns:
    None
    """

    # Read all the count columns from the files
    df_list = []
    for file in count_files:
        df = pd.read_csv(file, sep='\t', comment='#', index_col=0, usecols=[0, 6])
        df_list.append(df)

    # Read the gene annotation from the first file
    gene_annotation = pd.read_csv(count_files[0], sep='\t', comment='#', index_col=0, usecols=[0, 1, 2, 3, 4, 5])

    # If a gene metadata file is provided, merge it with the basic gene annotation
    merged_df = gene_annotation
    if genemeta_file:
        genemeta_df = pd.read_csv(genemeta_file, sep='\t', index_col=0)
        merged_df = merged_df.join(genemeta_df, how='left')

    # Concatenate all dataframes
    merged_df = merged_df.join(pd.concat(df_list, axis=1), how='left')

    # Remove ".markdup.sorted.bam" from all column names
    merged_df.columns = merged_df.columns.str.replace('.markdup.sorted.bam', '')

    # Write the merged dataframe to the output file
    merged_df.to_csv(output_file, sep='\t')

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Join multiple featureCounts matrices into one.")
    parser.add_argument("--genemeta", help="Gene metadata table to join.")
    parser.add_argument("--output", required=True, help="Output file for the merged matrix.")
    parser.add_argument("count_files", nargs='+', help="List of count matrix files to join.")

    args = parser.parse_args()

    join_featurecounts_matrix(args.count_files, args.output, args.genemeta)
