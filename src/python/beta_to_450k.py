#!/usr/bin/python3 -u

import argparse
import sys
import numpy as np
import os.path as op
import pandas as pd
from utils_wgbs import validate_single_file, validate_file_list, load_beta_data, \
                       beta2vec, IllegalArgumentError, eprint, ilmn2cpg_dict, \
                       add_multi_thread_args
from multiprocessing import Pool

# https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html


def single_beta(beta_path, indices, cov_thresh):
    return op.splitext(op.basename(beta_path))[0], \
           beta2vec(load_beta_data(beta_path)[indices - 1], min_cov=cov_thresh)


def read_reference(ref):
    # read Illumina-to-CpG_Index table:
    validate_single_file(ilmn2cpg_dict)
    df = pd.read_csv(ilmn2cpg_dict, sep='\t', header=None, names=['ilmn', 'cpg'])
    if ref is None:
        return df

    # validate and read reference file
    validate_single_file(ref)
    rf = pd.read_csv(ref, header=None, usecols=[0], names=['ilmn'])
    # remove first row if it's not a cg entry:
    if pd.isna(rf['ilmn'][0]) or not rf['ilmn'][0].startswith('cg'):
        rf = rf.iloc[1:, :]

    # merge reference file with map table
    mf = df.merge(rf, how='right', on='ilmn')

    # if there are sites that appear in the reference but not in the map table,
    # remove them and print a warning
    missing_sites = mf[mf['cpg'].isna()]
    if not missing_sites.empty:
        msg = 'WARNING: Skipping some unrecognized Illumina IDs \n'
        msg += f'(not found in the map table {ilmn2cpg_dict})\n'
        if not missing_sites['ilmn'].empty:
            eprint(missing_sites['ilmn'])
            eprint(list(missing_sites['ilmn']))
            msg += 'The missing sites: {}'.format(','.join(map(str, missing_sites['ilmn'])))
        eprint(msg)
    mf = mf[~mf['cpg'].isna()]

    mf['cpg'] = mf['cpg'].astype(int)
    return mf


def betas2csv(args):
    df = read_reference(args.ref)
    indices = np.array(df['cpg'])

    p = Pool(args.threads)
    params = [(b, indices, args.cov_thresh) for b in args.input_files]
    arr = p.starmap(single_beta, params)
    p.close()
    p.join()

    for name, beta in arr:
        df[name] = beta
    del df['cpg']

    df.to_csv(args.out_path, index=None, float_format='%.3f', na_rep='NA')


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+', help='one or more beta files')   # TODO: accept pat files (useful for small mixed files)
    parser.add_argument('-o', '--out_path', type=argparse.FileType('w'), default=sys.stdout,
                        help='Output path. [stdout]')
    parser.add_argument('-c', '--cov_thresh', type=int, default=3,
                        help='minimal coverage to include. sites with coverage < ' \
                        'cov_thresh are ignored [3]')
    parser.add_argument('--ref', help='a reference file with one column, ' \
                        'of Illumina IDs, optionally with a header line.')
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args


def main():
    """
    Convert beta file[s] to Illumina-450K format.
    Output: a csv file with ~480K rows, for the ~480K Illumina sites,
            and with columns corresponding to the beta files.
            all values are in range [0, 1], or NA.
            Only works for hg19.
    """
    args = parse_args()
    validate_file_list(args.input_files, '.beta')
    betas2csv(args)


if __name__ == '__main__':
    main()
