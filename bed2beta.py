#!/usr/bin/python3 -u

import argparse
import os.path as op
import sys
import pandas as pd
import numpy as np
from utils_wgbs import delete_or_skip, splitextgz, trim_to_uint8, validate_files_list, load_dict, eprint, load_dict_section
from genomic_region import GenomicRegion


def load_bed(bed_path, nrows, add1=False):
    df = pd.read_csv(bed_path, sep='\t', skiprows=0, header=None,
                     names=['chr', 'start', 'meth', 'total'],
                     usecols=[0, 1, 3, 4], nrows=nrows)
    nr_lines = df.shape[0]
    df.drop_duplicates(subset=['chr', 'start'], inplace=True)
    nr_lines_nodups = df.shape[0]
    if nr_lines != nr_lines_nodups:
        eprint('Warning: dropped {} duplicated lines'.format(nr_lines - nr_lines_nodups))
    if add1:
        df['start'] = df['start'] + 1
    return df


def bed2betas(args):

    # merge with the reference CpG bed file,
    # so the #lines in file will include all 28217448 sites (with NaN as 0)
    region = 'chr1:10469-876225' if args.debug else None
    try:
        rf = None       # Reference dictionary
        for bed in args.bed_paths:
            eprint('Converting {}...'.format(op.basename(bed)))
            # Check if bed should be skipped:
            outpath = op.join(args.outdir, splitextgz(op.basename(bed))[0]) + '.beta'
            if not delete_or_skip(outpath, args.force):
                continue

            # Load dict (at most once) and bed
            if rf is None:
                rf = load_dict_section(region, args.genome)
            df = load_bed(bed, nrows, args.add_one)

            # todo: multiprocess it, split by chromosome
            # merge dict with bed, then dump
            res = rf.merge(df, how='left', on=['chr', 'start']).fillna(0)
            trim_to_uint8(np.array(res[['meth', 'total']])).tofile(outpath)

    except pd.errors.ParserError as e:
        eprint('Invalid input file.\n{}'.format(e))
        return


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('bed_paths', nargs='+')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Overwrite existing files if existed')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--add_one', action='store_true',
                        help='Add 1 to start column, to match with CpG.bed.gz notation.')
                        #todo: infer that from bed.
    parser.add_argument('--outdir', '-o', default='.',
                        help='Output directory. Default is current directory [.]')
    parser.add_argument('--genome', help='Genome reference name. Default is hg19.',
                        default='hg19')
    args = parser.parse_args()
    return args


# def next_chrom(bed_path):   # todo: use or remove
#     try:
#         cur_chrom = ''
#         df = pd.DataFrame()
#         for chunk_df in pd.read_csv(bed_path, sep='\t', skiprows=0, header=None,
#                                     names=['chr', 'start', 'meth', 'total'],
#                                     usecols=[0, 1, 3, 4], chunksize=100000):
#             chrom = chunk_df['chr'][0]
#             cc = chunk_df[chunk_df['chr'] == cur_chrom]
#             nc = chunk_df[chunk_df['chr'] != cur_chrom]
#             if nc.shape[0]:
#                 cur_chrom = chrom
#             df = pd.concat([df, cc], reset_index=True)
#             if not df.empty():
#                 yield chunk_df
#             df = pd.concat([df, nc], reset_index=True)
#
#     except pd.errors.ParserError as e:
#         print('Invalid input file.\n{}'.format(e))
#         return


def main():
    """
    Convert bed[.gz] file[s] to beta file[s].
    bed file should be of the format (tab-separated):
    chr    start    end    #meth    #total
    """
    # todo: bed or bedGraph?
    args = parse_args()
    validate_files_list(args.bed_paths)
    bed2betas(args)


if __name__ == '__main__':
    main()

