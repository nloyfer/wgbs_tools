#!/usr/bin/python3 -u

import argparse
from utils_wgbs import add_GR_args, eprint, GenomeRefPaths, delete_or_skip, load_dict_section
from genomic_region import GenomicRegion
import pandas as pd
import numpy as np
from multiprocessing import Pool
import multiprocessing
import sys

COORDS_COLS = ['chr', 'start', 'end']


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    add_GR_args(parser, bed_file=True)
    parser.add_argument('--out_path', '-o', help='Output path for bed file [stdout]')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')
    parser.add_argument('-@', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='Number of threads to use (default: multiprocessing.cpu_count)')
    args = parser.parse_args()
    return args


def load_bed(bed_path, nrows):
    df = pd.read_csv(bed_path, sep='\t', header=None, nrows=nrows, comment='#')
    df.columns = COORDS_COLS + list(df.columns)[3:]
    # df.drop_duplicates(subset=['chr', 'start', 'end'], inplace=True)
    return df


def chr_thread(df, chrom, cf):
    # eprint(chrom)

    rf = load_dict_section(chrom).sort_values('start')

    # starts:
    s = pd.merge_asof(df.sort_values('start'), rf, by='chr', on='start', direction='forward')
    s = s.sort_values(by=['chr', 'start'])

    # ends:
    e = pd.merge_asof(df[['chr', 'end']].sort_values('end'), rf, by='chr', left_on='end', right_on='start',
                      direction='forward')
    e = e.sort_values(by=['chr', 'start'])
    e.loc[e['idx'].isna(), 'idx'] = cf[cf.chr == chrom]['size'].values[0]

    s['idx2'] = e['idx'].astype(int)
    s = s[COORDS_COLS + ['idx', 'idx2'] + list(s.columns)[3:-2]]
    return s


class AddCpGsToBed:
    def __init__(self, args):
        self.args = args
        self.out_path = args.out_path
        self.debug = args.debug
        if not delete_or_skip(self.out_path, self.args.force):
            return

        # load bed file:
        self.df = load_bed(args.bed_path, 100000 if self.debug else None)

        self.genome = GenomeRefPaths(args.genome)

        # load chromosomes sizes (in GpGs):
        self.cf = self.genome.get_chrom_cpg_size_table()
        self.cf['size'] = np.cumsum(self.cf['size'])
        self.proc_bed()

    def proc_bed(self):

        processes = []
        with Pool(self.args.threads) as p:
            chroms = [x for x in self.cf.chr if x in self.df.chr.unique()]
            for chrom in sorted(chroms):
                params = (self.df[self.df.chr == chrom], chrom, self.cf)
                processes.append(p.apply_async(chr_thread, params))
            p.close()
            p.join()

        r = pd.concat([pr.get() for pr in processes])
        if self.out_path is None:
            self.out_path = sys.stdout
        # r = self.reorder_to_original(r)
        r.to_csv(self.out_path, sep='\t', header=None, index=None, na_rep='NaN', mode='a')

    def reorder_to_original(self, res):
        return self.df[COORDS_COLS].merge(res, how='left', on=COORDS_COLS)


def convert_bed_file(args):
    """
    bed file should be of the format (tab-separated):
    Input:
    chr    start    end    [...]
    Output:
    chr    start    end    startCpG    endCpG   [...]
    """
    try:
        AddCpGsToBed(args)
    except pd.errors.ParserError as e:
        print('Invalid input file.\n{}'.format(e))
        return


def main():
    """
    Convert genomic region to CpG index range and vise versa
    """
    args = parse_args()

    if args.bed_path:
        convert_bed_file(args)
    else:
        print(GenomicRegion(args))


if __name__ == '__main__':
    main()
