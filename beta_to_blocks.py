#!/usr/bin/python3 -u

import argparse
import os
import numpy as np
import os.path as op
import pandas as pd
from utils_wgbs import validate_files_list
from multiprocessing import Pool
from os.path import splitext, basename
import sys
from utils_wgbs import load_beta_data, trim_to_uint8

# todo: supply this blocks file in wgbs_tools directory.
default_blocks_path = '/cs/cbio/netanel/blocks/outputs/nps20_genome.tsv.gz'


def apply_filter_wrapper(args, blocks_bins, finds, beta_path, df):

    try:
        # load beta file:
        data = load_beta_data(beta_path)

        # reduce to blocks:
        reduced_data = np.add.reduceat(data, blocks_bins)

        # filter data by min/max length and min/max #CpGs:
        reduced_data[:, 0][finds] = 0
        reduced_data[:, 1][finds] = 0

        # dump to file
        out_name = splitext(splitext(basename(args.blocks_file))[0])[0]
        out_name = splitext(basename(beta_path))[0] + '_' + out_name + '.bin'
        out_name = out_name.replace('_genome', '')
        out_name = op.join(args.out_dir, out_name)

        trim_to_uint8(reduced_data).tofile(out_name)
        print('saved to file:', out_name)

        if args.bedGraph:
            beta_vals = reduced_data[:, 0] / reduced_data[:, 1]
            # beta_vals[reduced_data[:, 1] == 0] = np.nan
            df['beta'] = beta_vals
            df.to_csv(out_name.replace('.bin', '.bedGraph'), sep='\t', index=None, header=None, na_rep=-1,
                      float_format='%.2f')

    except Exception as e:
        print('Failed with beta', beta_path)
        print('Exception:', e)


def main():
    """
    Collapse beta file to blocks binary file, of the same beta format
    """
    args = parse_args()
    files = args.input_files
    validate_files_list(files, '.beta')

    if not op.isfile(args.blocks_file):
        print('Invalid blocks file:', args.blocks_file)
        return

    names = ['chr', 'sloc', 'eloc', 'ssite', 'esite']
    df = pd.read_csv(args.blocks_file, sep='\t', usecols=[0, 1, 2, 3, 4], header=None, names=names)
    bplens = np.array(df['eloc'] - df['sloc']).flatten()
    cpglens = np.array(df['esite'] - df['ssite']).flatten()

    cpglen_indices = (args.max_sites < cpglens) | (cpglens < args.min_sites)
    bplen_indices = (args.max_bp < bplens) | (bplens < args.min_bp)
    filtered_indices = cpglen_indices | bplen_indices

    blocks_bins = np.array(df['ssite']).flatten() - 1

    print('starting Pool, with {} binary to collapse...'.format(len(files)))
    with Pool() as p:
        for beta_path in files:
            params = (args, blocks_bins, filtered_indices, beta_path, df[['chr', 'sloc', 'eloc']])
            p.apply_async(apply_filter_wrapper, params)
        p.close()
        p.join()

    # for beta_path in files:
    #     reduced_data = apply_filter_wrapper(beta_path, args.blocks, args.cov_thresh)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+', help='one or more beta files')
    parser.add_argument('-b', '--blocks_file', help='blocks path', default=default_blocks_path)
    parser.add_argument('-o', '--out_dir', help='output directory. Default is "."', default='.')
    parser.add_argument('--bedGraph', action='store_true', help='output a text file in addition to binary file')

    parser.add_argument('--min_sites', help='Minimum sites per block. Default is 0', type=int, default=0)
    parser.add_argument('--max_sites', help='Maximum sites per block. Default is inf', type=int, default=sys.maxsize)
    parser.add_argument('--min_bp', help='Minimal block length (bp). Default is 0', type=int, default=0)
    parser.add_argument('--max_bp', help='Maximal block length (bp). Default is inf', type=int, default=sys.maxsize)

    return parser.parse_args()


if __name__ == '__main__':
    main()
