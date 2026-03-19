#!/usr/bin/python3 -u

import argparse
import os.path as op
from utils_wgbs import load_beta_data, GenomeRefPaths, \
        trim_to_uint8, pretty_name, delete_or_skip


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('lbetas', nargs='+', help='one or more lbeta files')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if existed')
    parser.add_argument('-o', '--out_dir', help='Output directory for the beta file. [.]', default='.')
    parser.add_argument('--genome', help='Genome reference name.')
    args = parser.parse_args()
    return args

def main():
    """
    Generate a beta file from an lbeta file
    """
    args = parse_args()
    genome = GenomeRefPaths(args.genome)
    for lbeta in args.lbetas:
        out_beta = op.join(args.out_dir, pretty_name(lbeta)) + '.beta'
        if not delete_or_skip(out_beta, args.force):
            continue
        trim_to_uint8(load_beta_data(lbeta, genome=genome)).tofile(out_beta)


if __name__ == '__main__':
    main()
