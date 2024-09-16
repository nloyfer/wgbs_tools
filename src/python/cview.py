#!/usr/bin/python3 -u

import argparse
import subprocess
from utils_wgbs import MAX_PAT_LEN, pat_sampler, validate_single_file, \
    add_GR_args, cview_tool, collapse_pat_script, \
    cview_extend_blocks_script, validate_local_exe
from genomic_region import GenomicRegion
from beta_to_blocks import load_blocks_file


def subprocess_wrap_sigpipe(cmd):
    try:
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
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
        mpl = MAX_PAT_LEN
        if args.nanopore:
            mpl = 100000
        ms = max(1, s - mpl)
        cmd = f'tabix {pat} {gr.chrom}:{ms}-{e - 1} '

    view_flags = set_view_flags(args)
    cmd += f' | {cview_tool} --sites "{s}\t{e}" ' + view_flags
    cmd += add_subsample_cmd(args) # sub-sample reads
    if not gr.is_whole() and (('no_sort' not in args) or (not args.no_sort)):
        cmd += ' | sort -k2,2n -k3,3'
        if args.shuffle:
            cmd += 'R'
    cmd += f' | {collapse_pat_script} - '
    if get_cmd:
        return cmd
    if args.out_path is not None:
        cmd += f' > {args.out_path}'
    subprocess_wrap_sigpipe(cmd)


def add_subsample_cmd(args):
    if not hasattr(args, 'sub_sample') or args.sub_sample is None:
        return ''

    validate_local_exe(pat_sampler)
    ss = args.sub_sample
    rep = 1
    th = 0.25
    while ss > th:
        rep *= 2
        ss /= 2
    cmd = f' | {pat_sampler} {ss} {rep}'
    return cmd

def set_view_flags(args):
    view_flags = ''
    if args.strip:
        view_flags += ' --strip'
    if args.strict:
        view_flags += ' --strict'
    if args.min_len > 1:
        view_flags += f' --min_cpgs {args.min_len}'
    return view_flags


def view_bed(pat, args):
    # assume columns 4-5 of args.bed_file are startCpG, endCpG:
    bpath = args.bed_file

    # validate blocks file. If it's long, and starts with "chr1", use gunzip instead of tabix.
    df = load_blocks_file(bpath, nrows=1e6)
    if df.shape[0] == 1e6 and df.iloc[0, 0] in ('1', 'chr1'):
        tabix_cmd = f'gunzip -c {pat} '
    else:
        # extended blocks:
        tabix_cmd = 'gunzip -c' if bpath.endswith('.gz') else 'cat'
        tabix_cmd += f' {bpath} | {cview_extend_blocks_script} | tabix -R - {pat} '

    view_flags = set_view_flags(args)
    cmd = tabix_cmd + f' | {cview_tool} {view_flags} --blocks_path {bpath}'
    cmd += add_subsample_cmd(args) # sub-sample reads
    cmd += f' | sort -k2,2n -k3,3 | {collapse_pat_script} - '
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

def add_view_flags(parser, sub_sample=True, out_path=True, bed_file=True, long_reads=True):
    add_GR_args(parser, bed_file=bed_file)
    parser.add_argument('--strict', action='store_true',
                        help='pat: Truncate reads that start/end outside the given region. '
                             'Only relevant if "region", "sites" '
                             'or "bed_file" flags are given.')
    parser.add_argument('--strip', action='store_true',
                        help='pat: Remove trailing dots (from beginning/end of reads).')
    parser.add_argument('--min_len', type=int, default=1,
                        help='pat: Display only reads covering at least MIN_LEN CpG sites [1]')
    parser.add_argument('--shuffle', action='store_true',
                        help='pat: Shuffle reads order, while keeping the startCpG order '
                             '(sort -k2,2n -k3,3R)')
    parser.add_argument('--no_sort', action='store_true',
                        help='pat: Keep read order, as in the original pat file')
    if sub_sample:
        parser.add_argument('--sub_sample', type=float, #metavar='[0.0, 1.0]',
                            help='pat: subsample from reads. Only supported for pat')
    if out_path:
        parser.add_argument('-o', '--out_path', help='Output path. [stdout]')
    if long_reads:
        parser.add_argument('-np', '--nanopore', action='store_true',
                help='BETA VERSION: pull very long reads starting before the requested region')
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
    if (args.sub_sample is not None) and (args.sub_sample < 0):
        parser.error('[wt view] sub-sampling rate must be >= 0')
    validate_local_exe(cview_tool)
    cview(pat, args)



if __name__ == '__main__':
    main()
