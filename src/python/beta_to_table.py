#!/usr/bin/python3 -u

import argparse
import sys
import os.path as op
import warnings
from multiprocessing import Pool
import pandas as pd
import numpy as np
from dmb import load_gfile_helper, match_prefix_to_bin, load_uxm
from beta_to_blocks import collapse_process, load_blocks_file, is_block_file_nice
from utils_wgbs import validate_single_file, validate_file_list, eprint, \
    IllegalArgumentError, beta2vec, add_multi_thread_args, \
    drop_dup_keep_order


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'blocks', help='Blocks file with no header and with >= 5 columns')
    parser.add_argument(
        '--output', '-o', help='specify output path for the table [Default is stdout]')
    parser.add_argument('--groups_file', '-g',
                        help='groups csv file with at least 2 columns: name, group. beta files belong to the same group are averaged')
    parser.add_argument('--betas', nargs='+', help='beta files', required=True)
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('-c', '--min_cov', type=int, default=4, help='Minimal coverage to be considered. '
                        'blocks with less than MIN_COV site observations are considered as missing. [4]')
    parser.add_argument('--digits', type=int, default=2,
                        help='float percision (number of digits) [2]')
    parser.add_argument('--chunk_size', type=int, default=200000,
                        help='Number of blocks to load on each step [200000]')
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args


def groups_load_wrap(groups_file, betas):
    if groups_file is not None:
        validate_single_file(groups_file)
        validate_file_list(betas)
        gf = load_gfile_helper(groups_file)
    else:
        # otherwise generate dummy group file for all binary files in input_dir
        # first drop duplicated files, while keeping original order
        betas = drop_dup_keep_order(betas.copy())
        fnames = [op.splitext(op.basename(b))[0] for b in betas]
        gf = pd.DataFrame(columns=['fname'], data=fnames)
        gf['group'] = gf['fname']

    suff = '.beta'
    if betas[0].endswith('.lbeta'):
        suff = '.lbeta'
    gf['full_path'] = match_prefix_to_bin(gf['fname'], betas, suff)
    return gf


def cwrap(beta_path, blocks_df, is_nice, min_cov):
    if beta_path.endswith(('.beta', '.lbeta')):
        r = collapse_process(beta_path, blocks_df, is_nice)
        if r is None:
            return
        name = op.splitext(op.basename(beta_path))[0]
        return {name: beta2vec(r, min_cov)}
    else:
        eprint(f'[wt table] WARNING: {beta_path} is not a beta/lbeta file')

    return {op.basename(beta_path)[:-4]: load_uxm(beta_path, blocks_df, 'U', min_cov)}


def get_table(blocks_df, gf, min_cov, threads=8, verbose=False, group=True):
    is_nice, _ = is_block_file_nice(blocks_df)
    if verbose:
        eprint(f'[wt table] reducing to {blocks_df.shape[0]:,} blocks')
    betas = drop_dup_keep_order(gf['full_path'])
    p = Pool(threads)
    params = [(b, blocks_df, is_nice, min_cov) for b in betas]
    arr = p.starmap(cwrap, params)
    p.close()
    p.join()

    dicts = [d for d in arr if d is not None]
    dres = {k: v for d in dicts for k, v in d.items()}
    if not dres:
        fbetas = gf['fname'].tolist()
        eprint(f'[ wt table ] failed reducing {fbetas} to blocks\n{blocks_df}')
        raise IllegalArgumentError()

    if dres[list(dres.keys())[0]].size != blocks_df.shape[0]:
        eprint('[ wt table] beta2block returned wrong number of values')
        raise IllegalArgumentError()

    blocks_df.reset_index(drop=True, inplace=True)
    if not group:
        return pd.concat([blocks_df, pd.DataFrame(dres)[gf['fname'].tolist()]], axis=1)

    ugroups = drop_dup_keep_order(gf['group'])
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        empty_df = pd.DataFrame(index=blocks_df.index, columns=ugroups)
        blocks_df = pd.concat([blocks_df, empty_df], axis=1)
        for ugroup in ugroups:
            blocks_df[ugroup] = np.nanmean(
                np.concatenate([dres[k][None, :] for k in gf['fname'][gf['group'] == ugroup]]), axis=0).T
    return blocks_df


def betas2table(betas, blocks, groups_file, min_cov, threads=8, verbose=False):
    validate_single_file(blocks)
    gf = groups_load_wrap(groups_file, betas)
    blocks_df = load_blocks_file(blocks)
    return get_table(blocks_df, gf, min_cov, threads, verbose)


def dump(outpath, df, first=True, digits=3):
    if first:
        header = True
        mode = 'w'
    else:
        header = None
        mode = 'a'
    if outpath is None:
        outpath = sys.stdout
    df.to_csv(outpath, na_rep='NA',
              float_format=f"%.{digits}f",
              index=None, sep='\t',
              mode=mode, header=header)


def beta2table_generator(betas, blocks, groups_file, min_cov, threads, chunk_size=None, verbose=False):
    validate_single_file(blocks)
    gf = groups_load_wrap(groups_file, betas)
    blocks_df = load_blocks_file(blocks)
    if chunk_size is None:
        chunk_size = blocks_df.shape[0]
    for start in range(0, blocks_df.shape[0], chunk_size):
        subset_blocks = blocks_df.iloc[start:start + chunk_size].copy()
        yield get_table(subset_blocks, gf, min_cov, threads, verbose)


def main():
    """
    build a text table from beta files
    Optionally collapse samples with groups file
    """
    args = parse_args()
    chunks = beta2table_generator(args.betas,
                                  args.blocks,
                                  args.groups_file,
                                  args.min_cov,
                                  args.threads,
                                  args.chunk_size,
                                  args.verbose)
    first_chunk = True
    for chunk in chunks:
        dump(args.output, chunk, first_chunk, args.digits)
        first_chunk = False


if __name__ == '__main__':
    main()
