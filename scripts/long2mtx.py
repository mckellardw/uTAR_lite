import argparse
import pandas as pd
from scipy.sparse import coo_matrix
import os
import gzip


def parse_arguments():
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Convert a long format matrix to a Market Matrix (.mtx) or HDF5 (.h5) file"
    )
    parser.add_argument(
        "umitools_tsv",
        type=str,
        help="The input file containing the long format matrix",
    )
    parser.add_argument(
        "out_mat", type=str, help="The output file to save the matrix to"
    )
    parser.add_argument(
        "--output-format",
        type=str,
        choices=["mtx", "h5"],
        default="mtx",
        help="The format to save the output file in (default: mtx)",
    )
    return parser.parse_args()


def convert_to_mtx(umitools_tsv, out_mat, output_format):
    """
    Convert a long format matrix to a sparse matrix and save it in the specified format.

    Args:
        umitools_tsv (str): Path to the input file containing the long format matrix.
        out_mat (str): Path to the output file to save the matrix to.
        output_format (str): Format to save the output file in ('mtx' or 'h5').
    """
    if umitools_tsv.endswith(".gz"):
        # if input file is gzipped, use pandas read_csv with gzip compression
        long_df = pd.read_csv(
            umitools_tsv,
            compression="gzip",
            delimiter="\t",
            dtype={"gene": str, "cell": str, "count": int},
        )
    else:
        # if input file is not gzipped, use pandas read_csv without compression
        long_df = pd.read_csv(
            umitools_tsv, delimiter="\t", dtype={"gene": str, "cell": str, "count": int}
        )

    # Print summary of the loaded data
    print(f"Loaded data summary:")
    print(f"Number of unique genes: {long_df['gene'].nunique()}")
    print(f"Number of unique cells: {long_df['cell'].nunique()}")
    print(f"Total number of counts: {long_df['count'].sum()}")

    # save the 'gene' and 'cell' names as tsv.gz files
    gene_file = out_mat.replace(".mtx", "_genes.tsv.gz")
    cell_file = out_mat.replace(".mtx", "_cells.tsv.gz")
    long_df[["gene"]].drop_duplicates().to_csv(
        gene_file, sep="\t", index=False, header=False, compression="gzip"
    )
    long_df[["cell"]].drop_duplicates().to_csv(
        cell_file, sep="\t", index=False, header=False, compression="gzip"
    )

    # factorize the 'gene' and 'cell' names
    long_df["gene"] = pd.factorize(long_df["gene"])[0]
    long_df["cell"] = pd.factorize(long_df["cell"])[0]

    # convert the long format matrix to a sparse matrix
    sparse_mtx = coo_matrix((long_df["count"], (long_df["gene"], long_df["cell"])))

    if output_format == "mtx":
        import scipy.io as sio

        # save the sparse matrix as a .mtx file
        sio.mmwrite(out_mat, sparse_mtx)
    elif output_format == "h5":
        import h5py

        # save the sparse matrix as an .h5 file
        with h5py.File(out_mat, "w") as hf:
            hf.create_dataset(
                "utar_counts", data=sparse_mtx.data, compression="gzip", chunks=True
            )
            hf["utar_counts"].attrs["genes"] = sparse_mtx.row
            hf["utar_counts"].attrs["cells"] = sparse_mtx.col
    else:
        raise ValueError(f"Unsupported output format: {output_format}")


if __name__ == "__main__":
    args = parse_arguments()

    # check if the input file exists
    if not os.path.exists(args.umitools_tsv):
        raise FileNotFoundError(f"Input file '{args.umitools_tsv}' not found")

    convert_to_mtx(args.umitools_tsv, args.out_mat, args.output_format)
