#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_files_list, add_GR_args
from beta2bw import BetaToBigWig


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('beta_paths', nargs='+')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--outdir', '-o', default='.', help='Output directory. Default is current directory [.]')
    parser.add_argument('--remove_nan', action='store_true', help='If set, missing CpG sites are removed from the output.')
    add_GR_args(parser, bed_file = True)
    args = parser.parse_args()
    return args


def main():
    """
    Convert beta file[s] to bed file[s].
    """
    args = parse_args()
    validate_files_list(args.beta_paths, '.beta')
    b = BetaToBigWig(args)
    for beta in args.beta_paths:
        b.run_beta_to_bed(beta)


if __name__ == '__main__':
    main()

