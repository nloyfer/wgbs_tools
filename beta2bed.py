#!/usr/bin/python3 -u

import argparse
import os.path as op
import sys
import pandas as pd
from utils_wgbs import delete_or_skip, load_beta_data, load_dict, validate_files_list, eprint


def betas2beds(betas, outdir='.', genome='hg19', force=True, debug=False):

    # merge with the reference CpG bed file,
    # so the #lines in file will include all existing sites (with NaN as 0)
    nrows = 100000 if debug else None
    sites = (1, nrows + 1) if debug else None

    rf = None       # Reference dictionary
    for beta in betas:
        eprint('Converting {}...'.format(op.basename(beta)))
        # Check if beta should be skipped:
        outpath = op.join(outdir, op.splitext(op.basename(beta))[0]) + '.bed'   # todo: optional gzip?
        if not delete_or_skip(outpath, force):
            continue

        # Load dict (at most once) and beta
        if rf is None:
            rf = load_dict(nrows=nrows, genome_name=genome)
        barr = load_beta_data(beta, sites)

        # paste dict with beta, then dump
        rf['end'] = rf['start'] + 1
        rf['meth'] = barr[:, 0]
        rf['total'] = barr[:, 1]
        rf[rf['total'] > 0].to_csv(outpath, sep='\t', header=None, index=None)
        del rf['meth'], rf['total']



def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('beta_paths', nargs='+')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--outdir', '-o', default='.', help='Output directory. Default is current directory [.]')
    parser.add_argument('--genome', help='Genome reference name. Default is hg19.', default='hg19')

    args = parser.parse_args()
    return args


def main():
    """
    Convert beta file[s] to bed file[s].
    """
    args = parse_args()
    validate_files_list(args.beta_paths, '.beta')
    betas2beds(args.beta_paths, args.outdir, args.genome, args.force, args.debug)


if __name__ == '__main__':
    main()

