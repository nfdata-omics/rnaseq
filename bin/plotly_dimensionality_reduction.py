
from pathlib import Path
import argparse
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

def scatter_plot_with_metadata(input_table, xy_columns, xy_axis_labels):

    # read data with metadata
    df = pd.read_csv(input_table, index_col=0, sep="\t")

    # convert variables to categorical vars
    for mv in df.columns[2:]:
        df[mv] = df[mv].astype("category")

    # select metadata variables to include in the color list
    metavars = []
    for mv in df.columns[2:]:
        if len(df[mv].unique()) < 20:
            metavars.append(mv)

    # create plot and return it
    return nfdatautils.plt.scatter(
        df.reset_index(),
        xcolumn = xy_columns[0],
        ycolumn = xy_columns[1],
        samplenames = "samples",
        metavars = metavars,
        xaxis_label = xy_axis_labels[0],
        yaxis_lable = xy_axis_labels[1]
    )

def main():
    args = parse_args()

    # Check that the output folder exists
    if not args.output_folder.is_dir():
        raise FileNotFoundError(f"The specified output folder does not exist: {args.output_folder}")

    # define path of input and plot files
    data_path = args.output_folder / "downstream" / "dim_reduction" / "MDS_scores.txt"
    plot_path = args.output_folder / "downstream" / "dim_reduction" / "MDS_plot.html"

    # create plot and write it to html output file
    scatter_plot_with_metadata(data_path, xy_columns=("x", "y"), xy_axis_labels=("x_MDS", "y_MDS")).write_html(plot_path)

    # define path of input and plot files
    data_path = args.output_folder / "downstream" / "dim_reduction" / "PCA_scores.txt"
    plot_path = args.output_folder / "downstream" / "dim_reduction" / "PCA_plot.html"

    # create plot and write it to html output file
    scatter_plot_with_metadata(data_path, xy_columns=("PC1", "PC2"), xy_axis_labels=("PC1", "PC2")).write_html(plot_path)


if __name__ == "__main__":
    main()
