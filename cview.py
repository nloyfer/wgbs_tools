#!/usr/bin/python3 -u

import argparse
from utils_wgbs import load_beta_data2, MAX_PAT_LEN, pat_sampler, validate_single_file, \
    add_GR_args, IllegalArgumentError, BedFileWrap, load_dict_section, read_shell, eprint, add_multi_thread_args, \
    cview_tool, splitextgz, collapse_pat_script
from genomic_region import GenomicRegion
from convert import add_cpgs_to_bed, load_bed
from view import ViewPat
import subprocess
import uuid
import numpy as np
import sys
import os
import os.path as op
import pandas as pd

PAT_COLS = ('chr', 'start', 'pat', 'count')



###################
#                 #
#  Loading pat    #
#                 #
###################

def view_gr(pat, view_flags, args):
    gr = GenomicRegion(args)
    if gr.is_whole():
        s = 1
        e = 30000000   # todo: take last site from GenomeRef...
        cmd = f'gunzip -c {pat} '
    else:
        s, e = gr.sites
        s = max(1, s - MAX_PAT_LEN)
        cmd = f'tabix {pat} {gr.chrom}:{s}-{e - 1} '

    cmd = f"""/bin/bash -c 'cat <(echo "{s}\t{e}") <(echo "-1") <({cmd}) | {cview_tool}'"""
    cmd = cmd.rstrip("'") + f" {view_flags} ' "
    if args.sub_sample is not None:  # sub-sample reads
        cmd += f' | {pat_sampler} {args.sub_sample} '
    if not gr.is_whole():
        cmd += f' | sort -k2,2n -k3,3 '
    cmd += f' | {collapse_pat_script} - '
    # eprint(cmd)
    subprocess.check_call(cmd, shell=True)
    # cmd += f'wgbstools convert -L {bed} | cut -f4-5 | <(tabix -R - | {cview_tool}) '

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

def modified_bed_for_tabix(df, name):

    fake_df = df[['chr', 'startCpG', 'endCpG']].drop_duplicates()
    fake_df['startCpG'] = (fake_df['startCpG'].astype(int) - MAX_PAT_LEN).clip(lower=1)
    b = name + '.fake_tmp_df.tsv'
    fake_df.to_csv(b, sep='\t', header=None, index=None)
    t = b + '.tmp'
    cmd = f'sort -k1,1 -k2,2n {b} -o {t} && bedtools merge -d 1 -i {t} | sort -k2,2n > {b} && rm {t}'
    subprocess.check_call(cmd, shell=True)
    return b


def view_bed(pat, args):
    # df = add_cpgs_to_bed(args.bed_file, args.genome, drop_empty=True, threads=16)
    # assume columns 4-5 of args.bed_file are startCpG, endCpG:
    df = load_bed(args.bed_file).iloc[:, [0, 3, 4]]
    df.columns = ['chr', 'startCpG', 'endCpG']

    uniq_name = op.basename(splitextgz(pat)[0])
    uniq_name = f'{uniq_name}.tmp.{str(uuid.uuid4())[:8]}'
    uniq_name = op.join(args.tmp_dir, uniq_name)
    extended = modified_bed_for_tabix(df, uniq_name)
    view_flags = set_view_flags(args)
    # blocks_cmd = f'wgbstools convert -L {args.bed_file} | cut -f4-5 | sort -k1,1n '
    # assume columns 4-5 of args.bed_file are startCpG, endCpG:
    blocks_cmd = f'cut -f4-5 {args.bed_file} | sort -k1,1n '
    tabix_cmd = f'tabix -R {extended} {pat} '
    cmd = f"""/bin/bash -c 'cat <({blocks_cmd}) <(echo -1) <({tabix_cmd}) '"""
    cmd += f' | {cview_tool} {view_flags}'
    if args.sub_sample is not None:  # sub-sample reads
        cmd += f' | {pat_sampler} {args.sub_sample} '
    cmd += f' | sort -k2,2n -k3,3 | {collapse_pat_script} - '  # todo: implement the sort & collapsing in cview_tool
    # eprint(cmd)
    subprocess.check_call(cmd, shell=True)
    os.remove(extended)


def cview(pat, args):

    view_flags = set_view_flags(args)
    if args.bed_file:
        view_bed(pat, args)
    else:
        view_gr(pat, view_flags, args)

    return



##########################
#                        #
#         Main           #
#                        #
##########################


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat')
    add_GR_args(parser, bed_file=True)
    parser.add_argument('--tmp_dir', '-T', default='.',
            help='Temp directory for intermediate files [.]')
    parser.add_argument('--strict', action='store_true',
                        help='pat: Truncate reads that start/end outside the given region. '
                             'Only relevant if "region", "sites" '
                             'or "bed_file" flags are given.')
    parser.add_argument('--strip', action='store_true',
                        help='Remove trailing dots (from beginning/end of reads)')
    parser.add_argument('--min_len', type=int, default=1,
                        help='Display only reads covering at least MIN_LEN CpG sites [1]')
    parser.add_argument('--sub_sample', type=float, metavar='[0.0, 1.0]',
                        help='Subsample from reads')
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
