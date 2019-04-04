#!/usr/bin/python3 -u

import argparse
import sys
import numpy as np
import os.path as op
import pandas as pd
from utils_wgbs import validate_files_list, load_beta_data

ilmn2cpg_dict = '/cs/cbio/netanel/indexes/ilmn2CpG.tsv.gz'


def single_beta(beta_path, indices, cov_thresh):
    data = load_beta_data(beta_path)
    data = data[indices - 1]
    beta = np.divide(data[:, 0], data[:, 1], where=data[:, 1] >= cov_thresh)
    beta[data[:, 1] < cov_thresh] = np.nan
    return beta


def betas2csv(args):
    df = pd.read_table(ilmn2cpg_dict, header=None, names=['ilmn', 'cpg'])
    indices = np.array(df['cpg'])

    for beta_path in args.input_files:
        beta = single_beta(beta_path, indices, args.cov_thresh)
        name = op.splitext(op.basename(beta_path))[0]#[:20]
        df[name] = beta
    df = df.sort_values(by='ilmn')
    del df['cpg']
    df.to_csv(args.out_path, index=None, float_format='%.3f', na_rep='NaN')


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+', help='one or more beta files')
    parser.add_argument('-o', '--out_path', type=argparse.FileType('w'), default=sys.stdout,
                        help='Output path. [stdout]')
    parser.add_argument('-c', '--cov_thresh', type=int, default=1,
                        help='minimal coverage to include. sites with coverage < cov_thresh are ignored')
    args = parser.parse_args()
    return args


def main():
    """
    Convert beta file[s] to Illumina-450K format.
    Output: a csv file with ~480K rows, for the ~480K Illumina sites,
            and with columns corresponding to the beta files.
            all values are in range [0, 1], or NaN.
    """
    args = parse_args()
    validate_files_list(args.input_files, '.beta')
    betas2csv(args)


if __name__ == '__main__':
    main()
