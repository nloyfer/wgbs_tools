#!/usr/bin/python3 -u

import argparse
import sys
import os.path as op
from multiprocessing import Pool
import pandas as pd
import numpy as np
from utils_wgbs import validate_single_file, validate_file_list, load_beta_data, \
                       beta2vec, IllegalArgumentError, eprint, \
                       add_multi_thread_args, GenomeRefPaths, beta_sanity_check

# https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html


def single_beta(beta_path, indices, cov_thresh):
    name =  op.splitext(op.basename(beta_path))[0]
    values = beta2vec(load_beta_data(beta_path)[indices - 1], min_cov=cov_thresh)
    return name, values


def load_full_ref(args, genome):
    # read Illumina-to-CpG_Index table:
    df = pd.read_csv(validate_single_file(genome.ilmn2cpg_dict), sep='\t', header=None)

    # in old versions of ilmn2cpg_dict it only contains the 450K sites (no 450/850 column)
    # so ignore --EPIC flag, with a warning (not an error)
    if df.shape[1] < 3:
        if args.EPIC:
            msg = 'WARNING: current setup does not support translation' \
                        'to EPIC array.\n' \
                        'Try re-installing wgbstools and re-initializing current genome.' \
                        '\n--EPIC flag is ignored'
            eprint('[wt beta_450K] ' + msg)
        df['array'] = 450

    df.columns = ['ilmn', 'cpg', 'array']

    # filter to 450K sites (drop the 850K/EPIC ones) if not user ref file was provided:
    if (not args.EPIC) and (not args.ref):
        df = df[df['array'] == 450].reset_index(drop=True)

    df = df[['ilmn', 'cpg']]
    # If hg38, filter unmapped CpGs
    if genome.genome == 'hg38':
        df = df.dropna(how='any')
        df['cpg'] = df['cpg'].astype(int)
    return df


def read_reference(args):

    genome = GenomeRefPaths(args.genome)
    if not beta_sanity_check(args.input_files[0], genome):
        raise IllegalArgumentError('beta incompatible with genome')

    # load "full" reference - the one supplied with wgbstools
    df = load_full_ref(args, genome)
    if args.ref is None:
        return df

    # load user input reference file
    rf = pd.read_csv(validate_single_file(args.ref), header=None,
                     usecols=[0], names=['ilmn'])

    # remove first row if it's not a cg entry:
    # (This happens when the file has a header line)
    if pd.isna(rf['ilmn'][0]) or not rf['ilmn'][0].startswith('cg'):
        rf = rf.iloc[1:, :].reset_index(drop=True)

    # merge reference file with map table
    mf = df.merge(rf, how='right', on='ilmn')

    # if there are sites that appear in the reference but not in the map table,
    # remove them and print a warning
    missing_sites = mf[mf['cpg'].isna()]
    if not missing_sites.empty:
        msg = f'WARNING: Skipping {missing_sites.shape[0]} unrecognized Illumina IDs ' \
                f'(not found in the map table {genome.ilmn2cpg_dict})\n'
        if not missing_sites['ilmn'].empty:
            # eprint(missing_sites['ilmn'])
            # eprint(list(missing_sites['ilmn']))
            msg += 'Example for missing sites: {}'.format(','.join(map(str, missing_sites['ilmn'].head())))
            if not args.EPIC:
                msg += '\nconsider using --EPIC flag'
        eprint('[wt beta_450K] ' + msg)
    mf = mf[~mf['cpg'].isna()]

    mf['cpg'] = mf['cpg'].astype(int)
    return mf


def betas2csv(args):

    # set reference sites, as the intersection of the user input (--ref)
    # and the "full" reference, supplied by wgbstools (ilmn2cpg_dict)
    df = read_reference(args)
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
    parser.add_argument('--EPIC', action='store_true',
                        help='If no --ref is provided, output all 850K EPIC sites, instead of only 450K')
    add_multi_thread_args(parser)
    parser.add_argument('--genome', help='Genome reference name. Default is "default".', default='default')
    args = parser.parse_args()
    return args


def main():
    """
    Convert beta (or lbeta) file[s] to Illumina-450K (or EPIC) format.
    Output: a csv file with up to ~480K (or 850K) rows, for the Illumina array sites,
            and with columns corresponding to the beta files.
            all values are in range [0, 1], or NA.
            Only works for hg19 and hg38.
    """
    args = parse_args()
    if args.EPIC and args.ref:
        eprint('[wt beta_450K] WARNING: --ref is provided, therefore EPIC flag is ignored')
    validate_file_list(args.input_files)
    betas2csv(args)


if __name__ == '__main__':
    main()
