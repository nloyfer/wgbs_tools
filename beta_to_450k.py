#!/usr/bin/python3 -u

import argparse
import sys
import numpy as np
import os.path as op
import pandas as pd
from utils_wgbs import validate_files_list, load_beta_data, beta2vec, IllegalArgumentError, eprint
from multiprocessing import Pool, cpu_count

ilmn2cpg_dict = '/cs/cbio/netanel/indexes/ilmn2CpG.tsv.gz'  # todo: generate this and put in reference directory?


def single_beta(beta_path, indices, cov_thresh):
    return op.splitext(op.basename(beta_path))[0], \
           beta2vec(load_beta_data(beta_path)[indices - 1], min_cov=cov_thresh).astype(np.float16)


def read_reference(ref):
    # read Illumina-to-CpG_Index table:
    df = pd.read_csv(ilmn2cpg_dict, sep='\t', header=None, names=['ilmn', 'cpg'])
    if ref is None:
        return df

    # validate and read reference file
    if not op.isfile(ref):
        raise IllegalArgumentError('No such file: {}'.format(ref))
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
        msg += '(not found in the map table {})\n'.format(ilmn2cpg_dict)
        if not missing_sites['ilmn'].empty:
            print(missing_sites['ilmn'])
            print(list(missing_sites['ilmn']))
            msg += 'The missing sites: {}'.format(','.join(map(str, missing_sites['ilmn'])))
        eprint(msg)
    mf = mf[~mf['cpg'].isna()]

    mf['cpg'] = mf['cpg'].astype(int)
    return mf


def betas2csv(args):
    df = read_reference(args.ref)
    indices = np.array(df['cpg'])

    with Pool(cpu_count()) as p:
        resarr = [p.apply_async(single_beta, (b, indices, args.cov_thresh)) for b in args.input_files]
        p.close()
        p.join()

    for name, beta in [x.get() for x in resarr]:
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
                        help='minimal coverage to include. sites with coverage < cov_thresh are ignored [1]')
    parser.add_argument('--ref', help='a reference file with one column, of Illumina IDs, and no header.')
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
