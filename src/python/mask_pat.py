#!/usr/bin/python3 -u

import argparse
import subprocess
import os.path as op
from genomic_region import GenomicRegion
from pat2beta import pat2beta
from utils_wgbs import validate_single_file, mask_pat_tool, \
        delete_or_skip, main_script, collapse_pat_script, \
        add_GR_args, eprint, add_multi_thread_args


def mask_pat(pat_path, sites_to_hide, prefix, args):
    validate_single_file(pat_path)
    validate_single_file(sites_to_hide)

    pat_out = prefix + '.pat.gz'
    if not delete_or_skip(pat_out, args.force):
        return

    cmd = f'{main_script} cview {pat_path} '
    if args.bed_file:
        cmd += f' -L {args.bed_file}'
    else:
        gr = GenomicRegion(args)
        if not gr.is_whole():
            cmd += f' -r {gr.region_str}'
    cmd += f' | {mask_pat_tool} {args.sites_to_hide} '
    cmd += f' | {collapse_pat_script} -'
    cmd += f' | sort -k2,2n -k3,3'
    cmd += f' | bgzip -@ 4 > {pat_out}'
    cmd += f' && {main_script} index {pat_out}'
    subprocess.check_call(cmd, shell=True)
    if args.beta or args.lbeta:
        eprint(f'generate {"l" if args.lbeta else ""}beta file')
        pat2beta(pat_out, op.dirname(pat_out), args=args, force=True)
    return pat_out


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    add_GR_args(parser, bed_file=True)
    parser.add_argument('pat_path', help='pat file')
    add_multi_thread_args(parser)
    # mutually exclusive beta or lbeta:
    beta_or_lbeta = parser.add_mutually_exclusive_group(required=False)
    beta_or_lbeta.add_argument('--beta', action='store_true', help='create beta from output pat')
    beta_or_lbeta.add_argument('--lbeta', action='store_true', help='create lbeta from output pat')
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
