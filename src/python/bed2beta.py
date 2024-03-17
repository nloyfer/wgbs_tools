#!/usr/bin/python3 -u

import argparse
import os.path as op
import pandas as pd
import numpy as np
from utils_wgbs import delete_or_skip, splitextgz, trim_to_uint8, validate_file_list, \
                       eprint, load_dict_section


def load_bed(bed_path, nrows, add1=False):
    # check if there is a header:
    peek_df = pd.read_csv(bed_path, sep='\t', nrows=1, header=None)
    header = None if str(peek_df.iloc[0, 1]).isdigit() else 0
    df = pd.read_csv(bed_path, sep='\t', header=header,
                     names=['chr', 'start', 'meth', 'total'],
                     usecols=[0, 1, 3, 4], nrows=nrows)
    nr_lines = df.shape[0]
    df.drop_duplicates(subset=['chr', 'start'], inplace=True)
    nr_dup_lines = nr_lines - df.shape[0]
    if nr_dup_lines > 0:
        eprint(f'[wt bed] Warning: dropped {nr_dup_lines} duplicated lines')
    if add1:
        df['start'] = df['start'] + 1
    return df


def bed2betas(args):

    # merge with the reference CpG bed file,
    # so the #lines in file will include all 28217448 sites (with NaN as 0)
    region = 'chr1:10469-876225' if args.debug else None
    nrows = 10000 if args.debug else None
    try:
        rf = None       # Reference dictionary
        for bed in args.bed_paths:
            eprint(f'[wt bed] Converting {op.basename(bed)}...')
            # Check if bed should be skipped
            outpath = op.join(args.outdir, splitextgz(op.basename(bed))[0] + '.beta')
            if not delete_or_skip(outpath, args.force):
                continue

            # Load dict (at most once) and bed
            if rf is None:
                rf = load_dict_section(region, args.genome)
            df = load_bed(bed, nrows, args.add_one)

            # todo: implement in C++.
            # merge dict with bed, then dump
            res = rf.merge(df, how='left', on=['chr', 'start']).fillna(0)
            trim_to_uint8(np.array(res[['meth', 'total']])).tofile(outpath)

    except pd.errors.ParserError as e:
        eprint(f'[wt bed] Invalid input file.\n{e}')
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
    parser.add_argument('--genome', help='Genome reference name.')
    args = parser.parse_args()
    return args


def main():
    """
    Convert bed[.gz] file[s] to beta file[s].
    bed file should be of the format (tab-separated):
    chr    start    end    #meth    #total
    """
    args = parse_args()
    validate_file_list(args.bed_paths)
    bed2betas(args)


if __name__ == '__main__':
    main()
