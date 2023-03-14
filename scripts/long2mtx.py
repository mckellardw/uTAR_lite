import argparse
import pandas as pd
from scipy.sparse import coo_matrix, save_npz
import os
import gzip


def convert_to_mtx(input_file, output_file, output_format):
    if input_file.endswith('.gz'):
        # if input file is gzipped, use pandas read_csv with gzip compression
        long_df = pd.read_csv(
            input_file, 
            compression='gzip',
            delimiter="\t",
            dtype={'gene': str, 'cell': str, 'count': int}
        )
    else:
        # if input file is not gzipped, use pandas read_csv without compression
        long_df = pd.read_csv(
            input_file,
            delimiter="\t",
            dtype={'gene': str, 'cell': str, 'count': int}
        )

    # save the 'gene' and 'cell' names as tsv.gz files
    gene_file = output_file.replace('.mtx', '_genes.tsv.gz')
    cell_file = output_file.replace('.mtx', '_cells.tsv.gz')
    long_df[['gene']].drop_duplicates().to_csv(
        gene_file, 
        sep='\t', 
        index=False, 
        header=False,
        compression='gzip'
    )
    long_df[['cell']].drop_duplicates().to_csv(
        cell_file, 
        sep='\t', 
        index=False,
        header=False,
        compression='gzip'
    )

    # factorize the 'gene' and 'cell' names
    long_df['gene'] = pd.factorize(long_df['gene'])[0]
    long_df['cell'] = pd.factorize(long_df['cell'])[0]

    # convert the long format matrix to a sparse matrix
    sparse_mtx = coo_matrix((long_df['count'], (long_df['gene'], long_df['cell'])))

    if output_format == 'mtx':
        import scipy.io as sio
        # save the sparse matrix as a .mtx file
        sio.mmwrite(
            output_file, 
            sparse_mtx
        )
    elif output_format == 'h5':
        import h5py
        # save the sparse matrix as an .h5 file
        with h5py.File(output_file, 'w') as hf:
            hf.create_dataset(
                'utar_counts', 
                data=sparse_mtx.data, 
                compression='gzip',
                chunks=True# sparse_mtx.shape
            )
            hf['utar_counts'].attrs['genes'] = sparse_mtx.row
            hf['utar_counts'].attrs['cells'] = sparse_mtx.col
    else:
        raise ValueError(f"Unsupported output format: {output_format}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert a long format matrix to a Market Matrix (.mtx) or HDF5 (.h5) file'
    )
    parser.add_argument(
        'input_file', 
        type=str, 
        help='The input file containing the long format matrix'
    )
    parser.add_argument(
        'output_file', 
        type=str, 
        help='The output file to save the matrix to'
    )
    parser.add_argument(
        '--output-format', 
        type=str, 
        choices=['mtx', 'h5'], 
        default='mtx',
        help='The format to save the output file in (default: mtx)'
    )
    args = parser.parse_args()

    # check if the input file exists
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input file '{args.input_file}' not found")

    convert_to_mtx(args.input_file, args.output_file, args.output_format)
