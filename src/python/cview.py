#!/usr/bin/python3 -u

import argparse
from utils_wgbs import MAX_PAT_LEN, pat_sampler, validate_single_file, \
    add_GR_args, eprint, cview_tool, collapse_pat_script, \
    cview_extend_blocks_script, add_multi_thread_args
from genomic_region import GenomicRegion
from beta_to_blocks import load_blocks_file
import subprocess
import os.path as op


def subprocess_wrap_sigpipe(cmd):
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        if e.returncode != 141:   # e.g. if the output is piped to head
            raise e

###################
#                 #
#  Loading pat    #
#                 #
###################

def view_gr(pat, args, get_cmd=False):
    validate_single_file(pat, '.pat.gz')
    gr = GenomicRegion(args)
    if gr.is_whole():
        s = 1
        e = gr.genome.get_nr_sites() + 1
        cmd = f'gunzip -c {pat} '
    else:
        s, e = gr.sites
        ms = max(1, s - MAX_PAT_LEN)
        cmd = f'tabix {pat} {gr.chrom}:{ms}-{e - 1} '

    view_flags = set_view_flags(args)
    cmd = f"""/bin/bash -c 'cat <(echo "{s}\t{e}") <(echo "-1") <({cmd}) | {cview_tool}'"""
    cmd = cmd.rstrip("'") + f" {view_flags} ' "
    if hasattr(args, 'sub_sample') and args.sub_sample is not None:  # sub-sample reads
        cmd += f' | {pat_sampler} {args.sub_sample} '
    if not gr.is_whole():
        cmd += f' | sort -k2,2n -k3,3 '
    cmd += f' | {collapse_pat_script} - '
    if get_cmd:
        return cmd
    if args.out_path is not None:
        cmd += f' > {args.out_path}'
    subprocess_wrap_sigpipe(cmd)


def set_view_flags(args):
    view_flags = ''
    if args.strip:
        view_flags += ' --strip'
    if args.strict:
        view_flags += ' --strict'
    if args.min_len > 1:
        view_flags += f' --min_cpgs {args.min_len}'
    # if args.sub_sample:
        # view_flags += f' --sub_sample {args.sub_sample}'
    return view_flags


def view_bed(pat, args):
    # assume columns 4-5 of args.bed_file are startCpG, endCpG:
    bpath = args.bed_file

    # validate blocks file. If it's long, and starts with "chr1", use gunzip instead of tabix.
    tabix_cmd = ''
    df = load_blocks_file(bpath, nrows=1e6)
    if df.shape[0] == 1e6 and df.iloc[0, 0] in ('1', 'chr1'):
        tabix_cmd = f'gunzip -c {pat} '

    # extended blocks:
    cat_cmd = ('gunzip -c' if bpath.endswith('.gz') else 'cat') + f' {bpath}'
    if not tabix_cmd:
        tabix_cmd = cat_cmd + f' | {cview_extend_blocks_script} | tabix -R - {pat} '

    blocks_cmd = cat_cmd + f' | cut -f4-5 | sort -k1,1n '
    cmd = f'/bin/bash -c \"cat <({blocks_cmd}) <(echo -1) <({tabix_cmd}) \" '
    view_flags = set_view_flags(args)
    cmd += f' | {cview_tool} {view_flags}'
    if args.sub_sample is not None:  # sub-sample reads
        cmd += f' | {pat_sampler} {args.sub_sample} '
    cmd += f' | sort -k2,2n -k3,3 | {collapse_pat_script} - '
    # eprint(cmd)
    if args.out_path is not None:
        cmd += f' > {args.out_path}'
    subprocess_wrap_sigpipe(cmd)


def cview(pat, args):
    if args.bed_file:
        view_bed(pat, args)
    else:
        view_gr(pat, args)



##########################
#                        #
#         Main           #
#                        #
##########################

def add_view_flags(parser):
    add_GR_args(parser, bed_file=True)
    parser.add_argument('-o', '--out_path', help='Output path. [stdout]')
    parser.add_argument('--sub_sample', type=float, metavar='[0.0, 1.0]',
                        help='pat: subsample from reads. Only supported for pat')
    parser.add_argument('--strict', action='store_true',
                        help='pat: Truncate reads that start/end outside the given region. '
                             'Only relevant if "region", "sites" '
                             'or "bed_file" flags are given.')
    parser.add_argument('--strip', action='store_true',
                        help='pat: Remove trailing dots (from beginning/end of reads).')
    add_multi_thread_args(parser)
    parser.add_argument('--min_len', type=int, default=1,
                        help='pat: Display only reads covering at least MIN_LEN CpG sites [1]')
    return parser


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat')
    parser = add_view_flags(parser)
    return parser


def main():
    """ view pat file with the c++ engine """
    parser = parse_args()
    args = parser.parse_args()
    # validate input file
    pat = args.pat
    validate_single_file(pat)
    if args.sub_sample is not None and not 1 >= args.sub_sample >= 0:
        parser.error('[wt view] sub-sampling rate must be within [0.0, 1.0]')
    cview(pat, args)



if __name__ == '__main__':
    main()
