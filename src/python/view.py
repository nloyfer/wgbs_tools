#!/usr/bin/python3 -u

import argparse
import subprocess
import numpy as np
import sys
import os
from io import StringIO
import os.path as op
import pandas as pd
from multiprocessing import Pool
import multiprocessing
from utils_wgbs import load_beta_data2, MAX_PAT_LEN, pat_sampler, validate_single_file, \
    add_GR_args, IllegalArgumentError, BedFileWrap, read_shell, eprint, \
    add_multi_thread_args, catch_BrokenPipeError
from genomic_region import GenomicRegion

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
            cmd = f'gunzip -cd {self.pat_path} '
        else:
            start, end = self.gr.sites if sites is None else sites
            start = max(1, start - MAX_PAT_LEN)
            cmd = f'tabix {self.pat_path} '
            cmd += f'{self.gr.chrom}:{start}-{end - 1} '  # non-inclusive
        return cmd

    def strict_reads(self, df):
        # trim reads outside the gr
        start, end = self.gr.sites
        for idx, row in df.iterrows():
            rstart = row[1]
            pat = row[2]
            if rstart < start:
                df.loc[idx, 'pat'] = pat = pat[start - row['start']:]
                df.loc[idx, 'start'] = rstart = start
            if rstart + len(pat) > end:
                df.loc[idx, 'pat'] = pat[:end - df.loc[idx, 'start']]
        return df

    def sample_reads(self, df):
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

    def perform_view(self):
        # todo use for loop for large sections (or full file)
        df = read_shell(self.build_cmd(), names=get_pat_cols(self.pat_path))
        if df.empty:
            # eprint('empty')
            return df
        if self.gr.sites is not None:
            start, _ = self.gr.sites
            df = df[df['start'] + df['pat'].str.len() > start]

        if self.strict:                     # --strict
            df = self.strict_reads(df)
        if self.strip:                      # --strip
            df = self.strip_reads(df)
        if self.min_len > 1:                # --min_len
            df = df[df['pat'].str.len() >= self.min_len]
        if self.sub_sample:                 # --sub_sample
            self.sample_reads(df)
        df.reset_index(drop=True, inplace=True)
        return df


def get_pat_cols(pat_path):
    try:
        cols = list(PAT_COLS)
        peek = pd.read_csv(pat_path, sep='\t', nrows=1, header=None)
        # validate fields:
        chrom, site, pat, count = peek.values[0][:4]
        if not (str(site).isdigit() and str(count).isdigit() and set(pat) <= set('.CT')):
            eprint('[wt view] WARNING: Invalid first line in pat file:', peek.values)
        while len(peek.columns) > len(cols):
            cols += ['tag{}'.format(len(cols) - len(PAT_COLS) + 1)]
    except pd.errors.EmptyDataError as e:
        eprint('[wt view] WARNING: Empty pat file')
    return cols


def view_pat_mult_proc(input_file, strict, sub_sample,
        min_len, grs, i, step, strip, genome):
    reads = []
    cgrs = []
    for i in range(i, min(len(grs), i + step)):
        try:
            gr = GenomicRegion(region=grs[i], genome_name=genome)
            df = ViewPat(input_file, sys.stdout, gr, strict, sub_sample,
                    None, min_len, strip).perform_view()
            x = df.to_csv(sep='\t', index=None, header=None)
        except IllegalArgumentError as e:
            gr = grs[i] + ' - No CpGs'
            x = ''
        reads.append(x)
        cgrs.append(gr)
    return reads, cgrs


def is_bed_disjoint(b):
    cmd = f"""/bin/bash -c 'diff {b} <(bedtools intersect -a {b} -b {b} -wa)' > /dev/null """
    if subprocess.call(cmd, shell=True):
        eprint(f'[wt view] WARNING: bed file {b} regions are not disjoint.\n' \
                '                   Reads covering overlapping regions will be duplicated.\n' \
                '                   Use cview to avoid read duplication.')


def view_pat_bed_multiprocess(args, bed_wrapper):
    if not bed_wrapper:
        raise IllegalArgumentError('bed file is None')
    is_bed_disjoint(args.bed_file)

    full_regions_lst = list(bed_wrapper.fast_iter_regions())
    bigstep = 100
    for ch in range(0, len(full_regions_lst), bigstep):
        regions_lst = full_regions_lst[ch:ch + bigstep]

        n = len(regions_lst)
        # eprint(args.input_file, ch)
        step = max(1, n // args.threads)

        processes = []
        with Pool() as p:
            for i in range(0, n, step):
                params = (args.input_file, args.strict, args.sub_sample,
                        args.min_len, regions_lst, i, step,
                        args.strip, args.genome)
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


####################
#                  #
#  Loading beta    #
#                  #
####################

def get_beta_section(beta_path, gr):
    # load data
    data = load_beta_data2(beta_path, gr=gr.sites)
    # load loci
    cmd = f'tabix {gr.genome.revdict_path} {gr.chrom}:{gr.sites[0]}-{gr.sites[1]-1} | cut -f1-2'
    txt = subprocess.check_output(cmd, shell=True).decode()
    names = ['chr', 'start']
    df = pd.read_csv(StringIO(txt), sep='\t', header=None, names=names)
    df['end'] = df['start'] + 2
    df['meth'] = data[:, 0]
    df['total'] = data[:, 1]
    return df


def view_beta(beta_path, gr, opath, threads):
    """
    View beta file in given region/sites range
    :param beta_path: beta file path
    :param gr: a GenomicRegion object
    :param opath: output path (or stdout)
    """
    if not gr.is_whole():
        df = get_beta_section(beta_path, gr)
        df.to_csv(opath, sep='\t', index=None, header=None)
        return df
    # else, whole genome
    p = Pool(threads)
    chroms = list(gr.genome.get_chrom_cpg_size_table()['chr'])
    params = [(beta_path, GenomicRegion(region=c)) for c in chroms]
    arr = p.starmap(get_beta_section, params)
    p.close()
    p.join()
    # dump
    mode = 'w'
    for a in arr:
        a.to_csv(opath, sep='\t', index=None, mode=mode, header=None)
        mode = 'a'


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
    parser.add_argument('--print_region', action='store_true', help='pat: Prints region before reads')
    return parser


def main():
    """
    View the content of input file (pat/beta) as plain text.
    Possible filter by genomic region or sites range
    Output to stdout as default
    """
    parser = parse_args()
    args = parser.parse_args()

    if args.sub_sample is not None and not 1 >= args.sub_sample >= 0:
        parser.error('[wt view] sub-sampling rate must be within [0.0, 1.0]')

    # validate input file
    input_file = args.input_file
    validate_single_file(input_file)

    bed_wrapper = BedFileWrap(args.bed_file) if args.bed_file else None

    try:
        if op.splitext(input_file)[1] in ('.beta', '.lbeta', '.bin'):
            if bed_wrapper:
                eprint('Error: -L flag is not supported for beta files')  #TODO implement?
                exit(1)
            gr = GenomicRegion(args)
            view_beta(input_file, gr, args.out_path, args.threads)
        elif input_file.endswith('.pat.gz'):
            if bed_wrapper:
                view_pat_bed_multiprocess(args, bed_wrapper)
            else:
                from cview import view_gr
                view_gr(input_file, args)
        else:
            raise IllegalArgumentError('Unknown input format:', input_file)

    except BrokenPipeError:
        catch_BrokenPipeError()


if __name__ == '__main__':
    main()
