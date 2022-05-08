#!/usr/bin/python3 -u

import argparse
from utils_wgbs import delete_or_skip, validate_single_file, add_GR_args
from cview import subprocess_wrap_sigpipe
from genomic_region import GenomicRegion
from view import bview_build_cmd
import os


def beta2bed_build_cmd(beta_path, gr, bed_file, min_cov, mean, keep_na):
    cmd = bview_build_cmd(beta_path, gr, bed_file)
    if min_cov > 1:
        cmd += " | awk -v OFS='\\t' '{cov=$5;m=$4; if ($5<@COV) {cov=0;m=0};print $1,$2,$3,m,cov}'".replace('@COV', str(min_cov))
    if not keep_na:
        cmd += " | awk '$5>0'"
    if mean:
        cmd += " | awk '{r=-1; if($5>0){r=$4/$5}; printf \"%s\\t%d\\t%d\\t%.3g\\n\",$1,$2,$3,r}'"
    return cmd


def beta_to_bed(beta_path, gr, bed_file, min_cov, mean, keep_na, force, opath):
    if not delete_or_skip(opath, force):
        return

    cmd = beta2bed_build_cmd(beta_path, gr, bed_file, min_cov, mean, keep_na)
    if opath is not None:
        if opath.endswith('.gz'):
            cmd += ' | gzip -c '
        cmd += f' > {opath}'
    subprocess_wrap_sigpipe(cmd)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('beta_path')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Overwrite existing files if existed')
    parser.add_argument('--keep_na', action='store_true',
                        help='If set, missing CpG sites are not removed from the output')
    parser.add_argument('--mean', action='store_true',
                        help='Output a mean methylation value column instead of <meth,cov> columns')
    parser.add_argument('-c', '--min_cov', type=int, default=1,
                        help='Minimal coverage to consider when computing beta values.'
                             ' Default is 1 (include all observations). '
                             ' Sites with less than MIN_COV coverage are considered as missing.'
                             ' Only relevant if --mean is specified.')
    parser.add_argument('--outpath', '-o', default='/dev/stdout', help='Output path. [stdout]')
    add_GR_args(parser, bed_file=True)
    args = parser.parse_args()
    return args


def main():
    """
    Convert beta file to bed file.
    """
    args = parse_args()
    validate_single_file(args.beta_path) #, '.beta')
    gr = GenomicRegion(args)
    beta_to_bed(args.beta_path, gr, args.bed_file, args.min_cov, args.mean, args.keep_na, args.force, args.outpath)


if __name__ == '__main__':
    main()

