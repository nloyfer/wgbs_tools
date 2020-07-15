#!/usr/bin/python3 -u

import argparse
from utils_wgbs import load_beta_data2, MAX_PAT_LEN, MAX_READ_LEN, pat_sampler, validate_single_file, \
    add_GR_args, IllegalArgumentError, BedFileWrap, load_dict_section, read_shell, eprint
from genomic_region import GenomicRegion
import subprocess
import numpy as np
import sys
import os
import os.path as op
import pandas as pd
from multiprocessing import Pool
import multiprocessing

UNQ_COLS = ['chr', 'start', 'len', 'pat', 'count']
INF_UNQ_COLS = ['chr', 'idx', 'start', 'len', 'pat', 'count']
PAT_COLS = ('chr', 'start', 'pat', 'count')


###################
#                 #
#  Loading pat    #
#                 #
###################


class ViewPat:
    def __init__(self, pat_path, opath, gr, strict=False, sub_sample=None,
                 bed_wrapper=None, min_len=None, strip=False):
        self.pat_path = pat_path
        self.opath = opath
        self.min_len = min_len
        self.gr = gr
        self.strict = strict
        self.sub_sample = sub_sample
        self.bed_wrapper = bed_wrapper
        self.strip = strip

    def build_cmd(self, sites=None):
        """ Load a section from pat file using tabix """
        if not self.gr.chrom:  # entire pat file (no region filters)
            cmd = 'gunzip -cd {} '.format(self.pat_path)
        else:
            start, end = self.gr.sites if sites is None else sites
            start = max(1, start - MAX_PAT_LEN)
            cmd = 'tabix {} '.format(self.pat_path)
            cmd += '{}:{}-{} '.format(self.gr.chrom, start, end - 1)  # non-inclusive
        return cmd

    def trim_reads(self, df):
        # trim reads outside the gr
        if self.strict:
            start, end = self.gr.sites
            for idx, row in df.iterrows():
                rstart = row[1]
                pat = row[2]
                if rstart < start:
                    df.loc[idx, 'pat'] = pat = pat[start - row['start']:]
                    df.loc[idx, 'start'] = rstart = start
                if rstart + len(pat) > end:
                    df.loc[idx, 'pat'] = pat[:end - df.loc[idx, 'start']]
        if self.min_len > 1:
            df = df[df['pat'].str.len() >= self.min_len]
        return df

    def sample_reads(self, df):
        if self.sub_sample:  # sub-sample reads
            df['count'] = np.random.binomial(df['count'], self.sub_sample)
            df.drop(df[df['count'] == 0].index, inplace=True)
            df.reset_index(inplace=True, drop=True)

    def strip_reads(self, df):
        # Remove trailing dots from the right
        df['pat'] = df['pat'].str.rstrip('.')
        # Drop all dots reads
        df.drop(df[df.pat.str.len() == 0].index, inplace=True)

        # Remove trailing dots from the left
        def foo(row):
            pat = row[2]
            newpat = pat.lstrip('.')
            row[1] = int(row[1]) + len(pat) - len(newpat)
            row[2] = newpat
            return row

        cond = df['pat'].str.startswith('.')
        df.loc[cond] = df[cond].apply(foo, axis=1)
        df.sort_values(by=['start', 'pat'], inplace=True)
        return df

    def perform_view(self, dump=True):
        df = read_shell(self.build_cmd(),
                        names=get_pat_cols(self.pat_path))  # todo use for loop for large sections (or full file)
        if df.empty:
            eprint('empty')
            return df
        if self.gr.sites is not None:
            start, _ = self.gr.sites
            df = df[df['start'] + df['pat'].str.len() > start]

        df = self.trim_reads(df)
        if self.strip:
            df = self.strip_reads(df)
        self.sample_reads(df)
        df.reset_index(drop=True, inplace=True)
        if dump:
            df.to_csv(self.opath, sep='\t', index=None, header=None)
        return df

    def compose_awk_cmd(self):
        cmd = self.build_cmd()
        if self.strip:
            eprint('WARNING: strip flag not supported with awk_engine. Ignoring it')
        if self.gr.chrom:
            start, end = self.gr.sites
            cmd += ' | awk \'{if ($2 + length($3) > %s) {print;}}\' ' % start

            if self.strict:  # trim reads outside the gr
                cmd += ' | awk \'{(OFS="\t");s=$2; pat=$3;' \
                       ' if (s < %s) {s=%s;pat=substr(pat,%s-$2 + 1)}' \
                       ' if (s + length(pat) > %e) {pat=substr(pat, 0, %e-s)} ' \
                       ' print $1,s,pat,$4,$5,$6;}\' '.replace('%s', str(start)).replace('%e', str(end))
            if self.min_len > 1:
                cmd += ' | awk \'{(OFS="\t"); if (length($3) >= %m) {print $1,$2,$3,$4,$5,$6}}\' '.replace('%m', str(self.min_len))
        if self.sub_sample:  # sub-sample reads
            cmd += ' | {} {} '.format(pat_sampler, self.sub_sample)
        return cmd

    def perform_view_awk(self):
        subprocess.call(self.compose_awk_cmd(), shell=True, stdout=self.opath)

    def view_pat(self, awk=False):
        if awk:
            self.perform_view_awk()
        else:
            self.perform_view()


def get_pat_cols(pat_path):
    cols = list(PAT_COLS)
    peek = pd.read_csv(pat_path, sep='\t', nrows=1, header=None)
    # validate fields:
    chrom, site, pat, count = peek.values[0][:4]
    if not (str(site).isdigit() and str(count).isdigit() and set(pat) <= set('.CT')):
        eprint('WARNING: Invalid first line in pat file:', peek.values)
    while len(peek.columns) > len(cols):
        cols += ['tag{}'.format(len(cols) - len(PAT_COLS) + 1)]
    return cols


def view_pat_mult_proc(input_file, strict, sub_sample, min_len, grs, i, step, awk, strip):
    reads = []
    cgrs = []
    for i in range(i, min(len(grs), i + step)):
        gr = GenomicRegion(region=grs[i])
        if awk:
            try:
                cmd = ViewPat(input_file, sys.stdout, gr, strict, sub_sample, None, min_len, strip).compose_awk_cmd()
                x = subprocess.check_output(cmd, shell=True)
                if x:
                    x = x.decode()
            except IllegalArgumentError as e:
                gr = grs[i] + ' - No CpGs'
                x = ''
            # print('x', cmd, x)
        else:
            df = ViewPat(input_file, sys.stdout, gr, strict, sub_sample, None, min_len, strip).perform_view(dump=False)
            x = df.to_csv(sep='\t', index=None, header=None)
        reads.append(x)
        cgrs.append(gr)
    return reads, cgrs


def view_pat_bed_multiprocess(args, bed_wrapper):
    if not bed_wrapper:
        raise IllegalArgumentError('bed file is None')

    regions_lst = list(bed_wrapper.fast_iter_regions())
    n = len(regions_lst)
    step = max(1, n // args.threads)

    processes = []
    with Pool() as p:
        for i in range(0, n, step):
            params = (args.input_file, args.strict, args.sub_sample, args.min_len, regions_lst, i, step, args.awk_engine, args.strip)
            processes.append(p.apply_async(view_pat_mult_proc, params))
        p.close()
        p.join()
    # res = [sec.decode() for pr in processes for sec in pr.get()]
    for pr in processes:
        for reads, regions in zip(*pr.get()):
            if args.print_region:
                args.out_path.write(str(regions) + '\n')
            if not reads: # if the current region has no CpGs
                continue
            args.out_path.write(reads)


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
        return read_shell(self.build_cmd_unq(), names=UNQ_COLS)

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
        for chunk in pd.read_csv(self.unq_path, sep='\t', chunksize=chunksize, header=None, names=UNQ_COLS):
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
    data = load_beta_data2(beta_path, gr=gr.sites)
    np.savetxt(opath, data, fmt='%s', delimiter='\t')


##########################
#                        #
#         Main           #
#                        #
##########################


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_file')
    add_GR_args(parser, bed_file=True)
    parser.add_argument('-o', '--out_path', type=argparse.FileType('w'), default=sys.stdout,
                        help='Output path. [stdout]')
    parser.add_argument('--sub_sample', type=float, metavar='(0.0, 1.0)',
                        help='pat: subsample from reads. Only supported for pat')
    parser.add_argument('--strict', action='store_true',
                        help='pat: Truncate reads that start/end outside the given region. '
                             'Only relevant if "region", "sites" '
                             'or "bed_file" flags are given.')
    parser.add_argument('--strip', action='store_true',
                        help='pat: Remove trailing dots (from beginning/end of reads).\n'
                                   'Not supported with awk_engine')
    parser.add_argument('--inflate', action='store_true', help='unq: add CpG-Index column to the output')
    parser.add_argument('--awk_engine', action='store_true',
                        help='pat: use awk engine instead of python.\n'
                             'Its saves RAM when dealing with large regions.')
    # parser.add_argument('--multiprocess', '-@', type=int, default=16,
                        # help='pat: If bed file is specified, use multiple processors to read multiple.\n'
                             # 'regions in parallel. Default number of processors: 16.')
    parser.add_argument('-@', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='Number of threads to use (default: multiprocessing.cpu_count)')
    parser.add_argument('--min_len', type=int, default=1,
                        help='Pat: Display only reads covering at least MIN_LEN CpG sites [1]')
    parser.add_argument('--print_region', action='store_true', help='pat: Prints region before reads')
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
        eprint('sub-sampling rate must be within (0.0, 1.0)')
        return

    bed_wrapper = BedFileWrap(args.bed_file) if args.bed_file else None
    gr = GenomicRegion(args)

    try:
        if op.splitext(input_file)[1] in ('.beta', '.lbeta', '.bin'):
            if bed_wrapper:
                eprint('Warning: ingnoring -L flag')  #TODO implement?
            view_beta(input_file, gr, args.out_path)
        elif input_file.endswith('.pat.gz'):
            if bed_wrapper:
                view_pat_bed_multiprocess(args, bed_wrapper)
            else:
                vp = ViewPat(input_file, args.out_path, gr, args.strict, args.sub_sample, bed_wrapper, args.min_len, args.strip)
                vp.view_pat(args.awk_engine)
        elif input_file.endswith('.unq.gz'):
            grs = bed_wrapper.iter_grs() if bed_wrapper else [gr]
            for gr in grs:
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
