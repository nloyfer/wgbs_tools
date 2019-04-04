#!/usr/bin/python3 -u

import argparse
from utils_wgbs import load_beta_data, MAX_PAT_LEN, MAX_READ_LEN, pat_sampler, validate_single_file, \
    add_GR_args, IllegalArgumentError, BedFileWrap, load_dict_section, read_shell
from genomic_region import GenomicRegion
import subprocess
import numpy as np
import sys
import os
import pandas as pd

UNQ_COLS = ['chr', 'start', 'len', 'pat', 'count']
INF_UNQ_COLS = ['chr', 'idx', 'start', 'len', 'pat', 'count']
PAT_COLS = ['chr', 'start', 'pat', 'count']

###################
#                 #
#  Loading pat    #
#                 #
###################


class ViewPat:
    def __init__(self, pat_path, opath, gr, strict=False, sub_sample=None, bed_wrapper=None):
        self.pat_path= pat_path
        self.opath = opath
        self.gr = gr
        self.strict = strict
        self.sub_sample = sub_sample
        self.bed_wrapper = bed_wrapper

    def build_cmd(self):
        """ Load a section from pat file using tabix """
        if not self.gr.chrom:   # entire pat file (no region filters)
            cmd = 'gunzip -cd {} '.format(self.pat_path)
        else:
            start, end = self.gr.sites
            cmd = 'tabix {} '.format(self.pat_path)
            cmd += '{}:{}-{} '.format(self.gr.chrom, max(1, start - MAX_PAT_LEN), end - 1)  # non-inclusive
        return cmd

    def trim_reads(self, df):
        # trim reads outside the gr
        if self.strict:
            start, end = self.gr.sites
            for idx, row in df.iterrows():
                chrom, rstart, pat, count = row
                if rstart < start:
                    df.loc[idx, 'pat'] = pat = pat[start - row['start']:]
                    df.loc[idx, 'start'] = rstart = start
                if rstart + len(pat) > end:
                    df.loc[idx, 'pat'] = pat[:end - df.loc[idx, 'start']]

    def sample_reads(self, df):
        if self.sub_sample:  # sub-sample reads
            df['count'] = np.random.binomial(df['count'], self.sub_sample)
            df.drop(df[df['count'] == 0].index, inplace=True)

    def perform_view(self, dump=True):
        start, end = self.gr.sites
        df = read_shell(self.build_cmd(), header=None, names=PAT_COLS)  # todo use for loop for large sections (or full file)
        df = df[df['start'] + df['pat'].str.len() > start]

        self.trim_reads(df)
        self.sample_reads(df)
        df.reset_index(inplace=True, drop=True)
        if dump:
            df.to_csv(self.opath, sep='\t', index=None, header=None)
        return df

    def perform_view_awk(self):

        cmd = self.build_cmd()
        if self.gr.chrom:
            start, end = self.gr.sites
            cmd += ' | awk \'{if ($2 + length($3) > %s) {print;}}\' ' % start

            if self.strict:  # trim reads outside the gr
                cmd += ' | awk \'{(OFS="\t");s=$2; pat=$3;' \
                       ' if (s < %s) {s=%s;pat=substr(pat,%s-$2 + 1)}' \
                       ' if (s + length(pat) > %e) {pat=substr(pat, 0, %e-s)} ' \
                       ' print $1,s,pat,$4;}\' '.replace('%s', str(start)).replace('%e', str(end))

        if self.sub_sample: # sub-sample reads
            cmd += ' | {} {} '.format(pat_sampler, self.sub_sample)

        subprocess.call(cmd, shell=True, stdout=self.opath)

    def view_pat(self, awk=False):
        grs = self.bed_wrapper.iter_grs() if self.bed_wrapper else [self.gr]
        for gr in grs:
            self.gr = gr
            if awk:
                self.perform_view_awk()
            else:
                self.perform_view()

#################
#               #
#  Loading unq  #
#               #
#################


class ViewUnq:
    def __init__(self, unq_path, opath, gr, inflate):
        self.unq_path = unq_path
        self.opath = opath
        self.gr = gr
        self.inflate = inflate

    def build_cmd_unq(self):
        if not self.gr.chrom:
            return 'gunzip -cd {}'.format(self.unq_path)

        start, end = self.gr.bp_tuple
        cmd = 'tabix {} '.format(self.unq_path)
        cmd += '{}:{}-{} '.format(self.gr.chrom, max(1, int(start) - MAX_READ_LEN), end)
        cmd += ' | awk \'{if (($2 + $3) > %s) {print;}}\'' % start
        return cmd

    def load_sec_str(self):
        """ Read a section from unq file using tabix, direct it to opath """
        subprocess.call(self.build_cmd_unq(), shell=True, stdout=self.opath)

    def load_sec_df(self):
        """ Load a section from unq file using tabix into a DataFrame"""
        return read_shell(self.build_cmd_unq(), header=None, names=UNQ_COLS)

    def inflate_df(self, df):
        if df.empty:
            return
        # load relevant section from the dictionary:
        first_loc = df['start'].iloc[0]
        last_loc = df['start'].iloc[df.shape[0] - 1] + df['len'].iloc[df.shape[0] - 1]
        dict_region = '{}:{}-{}'.format(df['chr'].iloc[0], first_loc, last_loc)
        rf = load_dict_section(dict_region)

        # merge df with dictionary:
        res = pd.merge_asof(df, rf, by='chr', on='start', direction='forward')

        # dump
        res[INF_UNQ_COLS].to_csv(self.opath, sep='\t', index=None, header=None)

    def inflate_full_unq(self):
        chunksize = 10 ** 5
        for chunk in pd.read_table(self.unq_path, chunksize=chunksize, header=None, names=UNQ_COLS):
            for chrom in chunk['chr'].unique():
                self.inflate_df(chunk[chunk['chr'] == chrom])

    def view(self):
        if not self.inflate:
            return self.load_sec_str()

        if not self.gr.chrom:
            self.inflate_full_unq()
        else:
            self.inflate_df(self.load_sec_df())


####################
#                  #
#  Loading beta    #
#                  #
####################


def view_beta(beta_path, gr, opath):
    """
    View beta file in given region/sites range
    :param beta_path: beta file path
    :param gr: a GenomicRegion object
    :param opath: output path (or stdout)
    """
    data = load_beta_data(beta_path, gr.sites)
    np.savetxt(opath, data, fmt='%s', delimiter='\t')

##########################
#                        #
#         Main           #
#                        #
##########################


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_file')
    add_GR_args(parser)
    parser.add_argument('-o', '--out_path', type=argparse.FileType('w'), default=sys.stdout,
                        help='Output path. [stdout]')
    parser.add_argument('--sub_sample', type=float, metavar='(0.0, 1.0)',
                        help='subsample from reads. Only supported for pat')    # todo: support unq too
    parser.add_argument('-L', '--bed_file', help='Only output reads overlapping the input BED FILE')  # todo: support unq
    parser.add_argument('--strict', action='store_true', help='Truncate reads that start/end outside the given region. '
                                                              'Only relevant if "region" or "sites" flags are given.')
    parser.add_argument('--inflate', action='store_true', help='unq: add CpG-Index column to the output')
    parser.add_argument('--awk_engine', action='store_true', help='pat: use awk engine instead of python.\n'
                                                                  'Its saves RAM when dealing with large regions.')
    args = parser.parse_args()
    return args


def main():
    """
    View the content of input file (pat/unq/beta) as plain text.
    Possible filter by genomic region or sites range
    Output to stdout as default
    """
    args = parse_args()
    # validate input file
    input_file = args.input_file
    validate_single_file(input_file)

    if args.sub_sample is not None and not 1 > args.sub_sample > 0:
        print('sub-sampling rate must be within (0.0, 1.0)', file=sys.stderr)
        return

    if args.bed_file and (args.region or args.sites):
        print('-L, -s and -r are mutually exclusive', file=sys.stderr)
        return

    bed_wrapper = BedFileWrap(args.bed_file) if args.bed_file else None
    gr = GenomicRegion(args)

    try:
        if input_file.endswith('.beta'):
            view_beta(input_file, gr, args.out_path)
        elif input_file.endswith('.pat.gz'):
            vp = ViewPat(input_file, args.out_path, gr, args.strict, args.sub_sample, bed_wrapper)
            vp.view_pat(args.awk_engine)
        elif input_file.endswith('.unq.gz'):
            ViewUnq(input_file, args.out_path, gr, args.inflate).view()
        else:
            raise IllegalArgumentError('Unknown input format:', input_file)

    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE


if __name__ == '__main__':
    main()
