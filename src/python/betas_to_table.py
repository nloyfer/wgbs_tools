#!/usr/bin/python3 -u

import argparse
import subprocess
import numpy as np
import sys
import os.path as op
import pandas as pd
import math
import os
from multiprocessing import Pool
from dmb import load_gfile_helper, match_prefix_to_bin
from beta_to_blocks import collapse_process, load_blocks_file, is_block_file_nice
from utils_wgbs import validate_single_file, validate_file_list, eprint, \
                       IllegalArgumentError, beta2vec, add_multi_thread_args, drop_dup_keep_order




def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser = argparse.ArgumentParser()
    parser.add_argument('blocks', help='Blocks file with no header and with >= 5 columns')
    parser.add_argument('--output', '-o', help='specify output path for the table [Default is stdout]')
    parser.add_argument('--groups_file', '-g', help='csv file of groups')
    parser.add_argument('--betas', nargs='+', help='beta files', required=True)
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('-c', '--min_cov', type=int, default=4, help='Minimal coverage to be considered,'
                                                                 'In both groups. [4]')
    parser.add_argument('--digits', type=int, default=2, help='float percision (number of digits) [2]')
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args


def dump(outpath, df, digits, verbose):
    if outpath is None:
        df.to_csv(sys.stdout, na_rep='NA',
                  float_format=f"%.{digits}f",
                  index=None, sep='\t')
        return

    ixs = np.array_split(df.index, 100)
    header = True
    mode = 'w'
    eixs = enumerate(ixs)
    if verbose:
        eprint(f'[wt table] dumping table with shape {df.shape} to {outpath}')
        from tqdm import tqdm  # todo: drop if not installed
        eixs = tqdm(enumerate(ixs), total=len(ixs))
    for ix, subset in eixs:
        df.loc[subset].to_csv(outpath, na_rep='NA',
                              float_format=f"%.{digits}f",
                              index=None, sep='\t',
                              mode=mode, header=header)
        header = None
        mode = 'a'

def groups_load_wrap(groups_file, betas):
    if groups_file is not None:
        validate_single_file(groups_file)
        validate_file_list(betas)
        gf = load_gfile_helper(groups_file)
    else:
        # otherwise, generate dummy group file for all binary files in input_dir
        # first drop duplicated files, while keeping original order
        betas = drop_dup_keep_order(betas.copy())
        fnames = [op.splitext(op.basename(b))[0] for b in betas]
        gf = pd.DataFrame(columns=['fname'], data=fnames)
        gf['group'] = gf['fname']
    gf['full_path'] = match_prefix_to_bin(gf['fname'], betas, '.beta')
    return gf


def cwrap(beta_path, blocks_df, is_nice, verbose):
    if verbose:
        eprint('[wt table]', op.splitext(op.basename(beta_path))[0])
    return collapse_process(beta_path, blocks_df, is_nice)


def betas2table(betas, blocks, groups_file, min_cov, threads=8, verbose=False):
    validate_single_file(blocks)

    gf = groups_load_wrap(groups_file, betas)
    blocks_df = load_blocks_file(blocks)
    is_nice, _ = is_block_file_nice(blocks_df)
    if verbose:
        eprint(f'[wt table] reducing to {blocks_df.shape[0]:,} blocks')
    p = Pool(threads)
    # params = [(b, blocks_df, verbose) for b in sorted(gf['full_path'].unique())]
    params = [(b, blocks_df, is_nice, verbose) for b in drop_dup_keep_order(gf['full_path'])]
    arr = p.starmap(cwrap, params)
    p.close()
    p.join()

    dicts = [d for d in arr if d is not None]
    dres = {k: beta2vec(v, min_cov) for d in dicts for k, v in d.items()}
    # groups = sorted(gf['group'].unique())
    groups = drop_dup_keep_order(gf['group'])
    with np.warnings.catch_warnings():
        np.warnings.filterwarnings('ignore', r'Mean of empty slice')
        for group in groups:
            blocks_df[group] = np.nanmean(np.concatenate([dres[k][None, :] for k in gf['fname'][gf['group'] == group]]), axis=0).T
    return blocks_df


def main():
    """
    build a text table from beta files
    Optionally collapse samples with groups file
    """
    args = parse_args()

    df = betas2table(args.betas, args.blocks, args.groups_file,
                     args.min_cov, args.threads, args.verbose)
    dump(args.output, df, args.digits, args.verbose)


if __name__ == '__main__':
    main()

