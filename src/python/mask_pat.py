#!/usr/bin/python3 -u

import argparse
import subprocess
import os.path as op
import numpy as np
from genomic_region import GenomicRegion
from utils_wgbs import validate_single_file, mask_pat_tool, \
        delete_or_skip, main_script, \
        GenomeRefPaths, add_GR_args, eprint


def mask_pat(pat_path, sites_to_hide, prefix, args):
    validate_single_file(pat_path)
    validate_single_file(sites_to_hide)

    pat_out = prefix + '.pat.gz'
    if not delete_or_skip(pat_out, args.force):
        return

    # nr_sites = GenomeRefPaths(args.genome).get_nr_sites()
    cmd = f'{main_script} cview {pat_path} '
    if args.bed_file:
        cmd += f' -L {args.bed_file}'
    else:
        gr = GenomicRegion(args)
        if not gr.is_whole():
            cmd += f' -r {gr.region_str}'
    cmd += f'| {mask_pat_tool} {args.sites_to_hide} '
    cmd += f' | bgzip -@ 4 > {pat_out}'
    cmd += f' && {main_script} index {pat_out}'
    # eprint(cmd)
    subprocess.check_call(cmd, shell=True)
    return pat_out


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    add_GR_args(parser, bed_file=True)
    parser.add_argument('pat_path', help='pat file')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if existed')
    parser.add_argument('-p', '--prefix', help='Output prefix for the masked pat file', required=True)
    parser.add_argument('--sites_to_hide', '-b', help='bed file with sites / blocks to mask out.\n'
            'must be 5 column format (chr, start, end, startCpG, endCpG)', required=True)
    return parser.parse_args()


def main():
    """
    mask out sites from a pat file
    """
    args = parse_args()
    mask_pat(args.pat_path, args.sites_to_hide, args.prefix, args)


if __name__ == '__main__':
    main()
