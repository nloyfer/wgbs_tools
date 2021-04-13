#!/usr/bin/python3 -u

import argparse
import subprocess
import numpy as np
import sys
import os.path as op
import pandas as pd
import math
import os
from tqdm import tqdm
from dmb import load_gfile_helper, match_prefix_to_bin
from multiprocessing import Pool
from beta_to_blocks import collapse_process, load_blocks_file
from utils_wgbs import validate_single_file, validate_file_list, eprint, \
                       IllegalArgumentError, beta2vec, add_multi_thread_args




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
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args


def dump(outpath, df, verbose):
    if outpath is None:
        df.to_csv(sys.stdout, na_rep='NA',
                  float_format="%.2f",
                  index=None, sep='\t')
        return

    ixs = np.array_split(df.index, 100)
    header = True
    mode = 'w'
    eixs = enumerate(ixs)
    if verbose:
        eprint(f'dumping table to {outpath}')
        eixs = tqdm(enumerate(ixs), total=len(ixs))
    for ix, subset in eixs:
        df.loc[subset].to_csv(outpath, na_rep='NA',
                              float_format="%.2f",
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
        fnames = [op.splitext(op.basename(b))[0] for b in betas]
        gf = pd.DataFrame(columns=['fname'], data=fnames)
        gf['group'] = gf['fname']
    gf['full_path'] = match_prefix_to_bin(gf['fname'], betas)
    return gf

def cwrap(b, df, verbose):
    if verbose:
        eprint(op.splitext(op.basename(b))[0])
    return collapse_process(b, df)

def main():
    """
    build a text table from beta files
    Optionally collapse samples with groups file
    """
    args = parse_args()
    validate_single_file(args.blocks)

    gf = groups_load_wrap(args.groups_file, args.betas)
    df = load_blocks_file(args.blocks)
    p = Pool(args.threads)
    params = [(b, df, args.verbose) for b in sorted(gf['full_path'].unique())]
    arr = p.starmap(cwrap, params)
    p.close()
    p.join()

    dicts = [d for d in arr if d is not None]
    dres = {k: beta2vec(v, args.min_cov) for d in dicts for k, v in d.items()}
    groups = sorted(gf['group'].unique())
    for group in groups:
        df[group] = np.nanmean(np.concatenate([dres[k][None, :] for k in gf['fname'][gf['group'] == group]]), axis=0).T

    dump(args.output, df, args.verbose)


if __name__ == '__main__':
    main()

