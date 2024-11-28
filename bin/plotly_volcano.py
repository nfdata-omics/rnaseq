
from pathlib import Path
import argparse
import glob
import pandas as pd
import nfdatautils

def parse_args():
    parser = argparse.ArgumentParser(description="Parse output folder and metadata file paths.")

    # Argument for the output folder
    parser.add_argument(
        "--output_folder",
        type=Path,
        required=True,
        help="Path to the folder containing output files."
    )

    return parser.parse_args()

def volcano_plot_from_deseq2_toptable(input_table):

    # read data with metadata
    df = pd.read_csv(input_table, index_col=0, sep="\t")

    # create plot and return it
    return nfdatautils.plt.volcano(
        df.reset_index(),
        lfc="log2FoldChange",
        pv="padj",
        axxlabel=r'log_2(FC)',
        axylabel=r'-log_10(pvalue)',
        dim=(10, 6),
        dotsize=10,
        valpha=None,
        geneid="gene_name"
    )

def main():
    args = parse_args()

    # Check that the output folder exists
    if not args.output_folder.is_dir():
        raise FileNotFoundError(f"The specified output folder does not exist: {args.output_folder}")

    # get the list of all dea tables present in the output folder
    deseq2_table_list = glob.glob(str(args.output_folder / "downstream/dea_*/deseq2_toptable.*.txt"))

    for data_path in deseq2_table_list:

        # define path of plot file
        plot_path = str(data_path).removesuffix(".txt") + ".html"

        # create plot and write it to html output file
        volcano_plot_from_deseq2_toptable(data_path).write_html(plot_path)


if __name__ == "__main__":
    main()

