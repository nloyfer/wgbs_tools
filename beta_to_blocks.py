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
from utils_wgbs import load_beta_data, trim_to_uint8, default_blocks_path, eprint


def get_bins(df):
    end = 28217449  # todo: read this from a reference file
    arr = np.unique(np.concatenate([[1], df['ssite'], df['esite'], [end]]))
    arr.sort()
    isin = np.isin(arr, np.concatenate([df['ssite'], [df['esite'][df.shape[0] - 1]]]))
    return arr - 1, isin


def apply_filter_wrapper(args, blocks_bins, finds, beta_path, df):
    try:
        # load beta file:
        data = load_beta_data(beta_path)

        # reduce to blocks:
        blocks_bins[-1] -= 1
        reduced_data = np.add.reduceat(data, blocks_bins)[finds][:-1]

        # dump to file
        out_name = splitext(splitext(basename(args.blocks_file))[0])[0]
        out_name = splitext(basename(beta_path))[0] + '_' + out_name + '.bin'
        out_name = out_name.replace('_genome', '')
        out_name = op.join(args.out_dir, out_name)

        trim_to_uint8(reduced_data).tofile(out_name)
        print('saved to file:', out_name)

        if args.bedGraph:
            with np.errstate(divide='ignore', invalid='ignore'):
                beta_vals = reduced_data[:, 0] / reduced_data[:, 1]
                print(beta_vals.shape, df.shape)
            # beta_vals[reduced_data[:, 1] == 0] = np.nan
            df['beta'] = beta_vals
            df.to_csv(out_name.replace('.bin', '.bedGraph'), sep='\t',
                      index=None, header=None, na_rep=-1,
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

    nr_removed = df[df.ssite == df.esite].shape[0]
    if nr_removed:
        eprint('removed {} regions with no CpGs'.format(nr_removed))

    if args.debug:
        print(df[df.ssite == df.esite])

    df = df[df.ssite < df.esite]
    blocks_bins, filtered_indices = get_bins(df)

    with Pool() as p:
        for beta_path in files:
            params = (args, blocks_bins,
                      filtered_indices, beta_path, df[['chr', 'sloc', 'eloc']])
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
    parser.add_argument('--debug', '-d', action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    main()
