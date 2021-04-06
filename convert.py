#!/usr/bin/python3 -u

import argparse
from utils_wgbs import add_GR_args, eprint, GenomeRefPaths, delete_or_skip, \
                       load_dict_section, add_multi_thread_args
from genomic_region import GenomicRegion
import pandas as pd
import numpy as np
from multiprocessing import Pool
import sys
import subprocess
from io import StringIO

COORDS_COLS = ['chr', 'start', 'end']


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    add_GR_args(parser, bed_file=True, no_anno=True)
    parser.add_argument('--out_path', '-o', help='Output path for bed file [stdout]')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--drop_empty', action='store_true', help='Drop empty regions (without CpGs)')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args


def load_bed(bed_path, nrows=None):
    df = pd.read_csv(bed_path, sep='\t', header=None, nrows=nrows, comment='#')
    df.columns = COORDS_COLS + list(df.columns)[3:]
    return df


def chr_thread(df, chrom, cf, genome):
    # eprint(chrom)

    # we don't need duplicates here (they'll be back later)
    df.drop_duplicates(subset=COORDS_COLS, inplace=True)
    rf = load_dict_section(chrom, genome).sort_values('start')

    # starts:
    s = pd.merge_asof(df.sort_values('start'), rf, by='chr', on='start', direction='forward')
    s = s.sort_values(by=['chr', 'start'])

    # ends:
    e = pd.merge_asof(df[['chr', 'end']].sort_values('end'), rf, by='chr', left_on='end', right_on='start',
                      direction='forward')
    e = e.sort_values(by=['chr', 'start'])
    e.loc[e['idx'].isna(), 'idx'] = cf[cf.chr == chrom]['size'].values[0] + 1

    # astype 'Int64' and not Int to avoid implicit casting to float in case there are NaNs
    s['idx2'] = e['idx'].astype('Int64')
    s['idx'] = s['idx'].astype('Int64')
    s = s[COORDS_COLS + ['idx', 'idx2'] + list(s.columns)[3:-2]]

    # drop regions without CpGs (pd.merge left them with a fictive 0 range, e.g 25-25)
    s.rename(columns={'idx': 'startCpG', 'idx2': 'endCpG'}, inplace=True)
    s.dropna(inplace=True, subset=['startCpG', 'endCpG'])
    s = s[s['endCpG'] - s['startCpG'] > 0]
    return s


def add_cpgs_to_bed(bed_file, genome, drop_empty, threads, add_anno):
    # load bed file:
    df = load_bed(bed_file)

    # load chromosomes sizes (in GpGs):
    cf = GenomeRefPaths(genome).get_chrom_cpg_size_table()
    cf['size'] = np.cumsum(cf['size'])

    chroms = sorted(set(cf.chr) & set(df.chr))
    params = [(df[df.chr == chrom].copy(), chrom, cf, genome) for chrom in chroms]
    p = Pool(threads)
    arr = p.starmap(chr_thread, params)
    p.close()
    p.join()

    # concat chromosomes
    r = pd.concat(arr)
    # merge with original table, to keep the order and the empty regions
    r = df.merge(r[COORDS_COLS + ['startCpG', 'endCpG']], how='left', on=COORDS_COLS)
    r = r[COORDS_COLS + ['startCpG', 'endCpG'] + list(r.columns)[3:-2]]

    # add annotations:
    if add_anno:
        annodf = get_anno(bed_file, genome)
        if annodf is not None:
            r = r.merge(annodf, how='left', on=COORDS_COLS)

    # drop regions w/o CpGs
    if drop_empty:
        r.dropna(inplace=True, subset=['startCpG', 'endCpG'])

    return r


def get_anno(bed_path, genome):
    anno_path = GenomeRefPaths(genome).annotations
    if anno_path is None:
        return
    read = 'gunzip -cd' if bed_path.endswith('.gz') else 'cat'
    cmd = f'{read} {bed_path} | cut -f1-3 | '
    cmd += f'bedtools intersect -wao -a - -b {anno_path} '
    cmd += '| bedtools merge -c 7,8 -o distinct,distinct '
    txt = subprocess.check_output(cmd, shell=True).decode()
    names = COORDS_COLS + ['type', 'gene']
    annodf = pd.read_csv(StringIO(txt), sep='\t', header=None, names=names)
    return annodf


def convert_bed_file(args):
    """
    bed file should be of the format (tab-separated):
    Input:
    chr    start    end    [...]
    Output:
    chr    start    end    startCpG    endCpG   [...]
    """
    try:
        out_path = args.out_path
        if not delete_or_skip(out_path, args.force):
            return
        if out_path is None:
            out_path = sys.stdout
        # add CpG columns
        r = add_cpgs_to_bed(bed_file=args.bed_file,
                genome=args.genome,
                drop_empty=args.drop_empty,
                threads=args.threads,
                add_anno=not args.no_anno)
        r.to_csv(out_path, sep='\t', header=None, index=None, na_rep='NaN')

    except pd.errors.ParserError as e:
        print(f'Invalid input file.\n{e}')
        return


def main():
    """
    Convert genomic region to CpG index range and vise versa
    """
    args = parse_args()

    if args.bed_file:
        convert_bed_file(args)
    else:
        print(GenomicRegion(args))


if __name__ == '__main__':
    main()
