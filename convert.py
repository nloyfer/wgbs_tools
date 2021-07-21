#!/usr/bin/python3 -u

import argparse
from utils_wgbs import add_GR_args, eprint, GenomeRefPaths, delete_or_skip, \
                       load_dict_section, add_multi_thread_args, IllegalArgumentError
from genomic_region import GenomicRegion
import pandas as pd
import numpy as np
from multiprocessing import Pool
import sys
import subprocess
from io import StringIO

COORDS_COLS3 = ['chr', 'start', 'end']
COORDS_COLS5 = COORDS_COLS3 + ['startCpG', 'endCpG']


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    region_or_sites = add_GR_args(parser, bed_file=True, no_anno=True)
    region_or_sites.add_argument('--site_file',
                        help='text file with a single CpG indexes column,' \
                             ' or <startCpG, endCpG> columns.\n' \
                             'if "-" is passed, the file is read from stdin.')
    parser.add_argument('--out_path', '-o', help='Output path for bed file [stdout]')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-p', '--parsable', action='store_true',
                        help='Output a parsing friendly format (only work with -r/-s flags)')
    parser.add_argument('--drop_empty', action='store_true',
                        help='Drop empty regions (without CpGs)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Overwrite existing files if existed')
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args

####################################################################
#                                                                  #
#                      Add CpGs to Bed                             #
#                                                                  #
####################################################################

def convert_bed_file(args):
    """
    bed file should be of the format (tab-separated):
    Input:
    chr    start    end    [...]
    Output:
    chr    start    end    startCpG    endCpG   [...]
    """
    out_path = sys.stdout if args.out_path is None else args.outpath
    if not delete_or_skip(out_path, args.force):
        return
    # add CpG columns
    bed_file = args.bed_file
    if bed_file == '-':
        bed_file = sys.stdin
    r = add_cpgs_to_bed(bed_file=bed_file,
            genome=args.genome,
            drop_empty=args.drop_empty,
            threads=args.threads,
            add_anno=not args.no_anno)
    r.to_csv(out_path, sep='\t', header=None, index=None, na_rep='NA')


def load_bed(bed_path, nrows=None):
    df = pd.read_csv(bed_path, sep='\t', header=None, nrows=nrows, comment='#')
    df.columns = COORDS_COLS3 + list(df.columns)[3:]
    return df


def slow_conversion(df, genome):
    df = df.iloc[:, :3]
    startCpGs = []
    endCpGs = []
    for ind, row in df.iterrows():
        try:
            sites = GenomicRegion(region='{}:{}-{}'.format(*row), genome_name=genome).sites
        except IllegalArgumentError as e:
            sites = (np.nan, np.nan)
        startCpGs.append(sites[0])
        endCpGs.append(sites[1])
    df['startCpG'] = pd.Series(startCpGs, dtype='Int64').values
    df['endCpG'] = pd.Series(endCpGs, dtype='Int64').values
    return df

def chr_thread(df, cf, genome):
    chrom = df['chr'].values[0]

    # we don't need duplicates here (they'll be back later)
    df.drop_duplicates(subset=COORDS_COLS3, inplace=True)

    # there must not be overlaps between blocks
    # if so, fall back to region-by-region conversion
    tf = df.sort_values(by='start')
    if not (tf['start'][1:].values - tf['end'][:tf.shape[0] - 1].values  >= 0).all():
        if df.shape[0] > 30:
            eprint(f'[wt convert] [{chrom}] WARNING: ' \
                    'Found overlaps in the input bed file. Conversion may be slow.')
        return slow_conversion(df, genome)

    rf = load_dict_section(chrom, genome).sort_values('start')

    # starts:
    s = pd.merge_asof(df.sort_values('start'), rf, by='chr', on='start', direction='forward')
    s = s.sort_values(by=['chr', 'start'])

    # ends:
    e = pd.merge_asof(df[['chr', 'end']].sort_values('end'), rf,
            by='chr', left_on='end', right_on='start', direction='forward')
    e = e.sort_values(by=['chr', 'start'])
    e.loc[e['idx'].isna(), 'idx'] = cf[cf.chr == chrom]['size'].values[0] + 1
    e.loc[(e.end == e.start).values, 'idx'] =  e.loc[(e.end == e.start).values, 'idx'] + 1

    # astype 'Int64' and not Int to avoid implicit casting to float in case there are NaNs
    s['idx2'] = e['idx'].astype('Int64')
    s['idx'] = s['idx'].astype('Int64')
    s = s[COORDS_COLS3 + ['idx', 'idx2'] + list(s.columns)[3:-2]]

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
    params = [(df[df.chr == chrom].copy(), cf, genome) for chrom in chroms]
    p = Pool(threads)
    arr = p.starmap(chr_thread, params)
    p.close()
    p.join()

    # concat chromosomes
    r = pd.concat(arr)
    # merge with original table, to keep the order and the empty regions
    r = df.merge(r[COORDS_COLS5], how='left', on=COORDS_COLS3)
    r = r[COORDS_COLS5 + list(r.columns)[3:-2]]

    # add annotations:
    if add_anno:
        annodf = get_anno(r, bed_file, genome)
        if annodf is not None:
            r = r.merge(annodf, how='left', on=COORDS_COLS3)

    # drop regions w/o CpGs
    if drop_empty:
        r.dropna(inplace=True, subset=['startCpG', 'endCpG'])

    return r


####################################################################
#                                                                  #
#                      Add Loci to CpGs                            #
#                                                                  #
####################################################################


def convert_site_file(args):
    """
    site file should be of the format (tab-separated):
    Input:
    startCpG    [endCpG]
    Output:
    chr    start    end    startCpG    [endCpG]
    """
    out_path = sys.stdout if args.out_path is None else args.outpath
    if not delete_or_skip(out_path, args.force):
        return
    # add loci columns
    r = add_bed_to_cpgs(cpgs=load_site_file(args.site_file),
                        genome=args.genome,
                        threads=args.threads)
                        # add_anno=not args.no_anno)
    r.to_csv(out_path, sep='\t', header=None, index=None, na_rep='NA')


def load_site_file(site_file):
    # load site file
    if site_file == '-':
        site_file = sys.stdin
    cpgs = pd.read_csv(site_file, sep='\t', header=None)
    nr_cols = cpgs.shape[1]
    if nr_cols == 1:
        cpgs.columns = ['startCpG']
        cpgs['endCpG'] = cpgs['startCpG'] + 1
    elif nr_cols == 2:
        cpgs.columns = ['startCpG', 'endCpG']
    else:
        raise IllegalArgumentError('Invalid site file. required columns: startCpG    [endCpG]')
    if not (cpgs['startCpG'] < cpgs['endCpG']).all():
        raise IllegalArgumentError('Invalid site file. endCpG must be greater thatn startCpG')
    return cpgs

def chr_thread_cpg(df, chrom, genome):
    # load reference dict:
    dict_df = load_dict_section(chrom, genome)
    df.drop_duplicates(inplace=True)

    # merge startCpG
    df = df.merge(dict_df, left_on='startCpG', right_on='idx')

    # merge endCpG
    dict_df = dict_df.rename(columns={'idx': 'endCpG', 'start': 'end'})
    df['endCpG'] = df['endCpG'] - 1
    df = df.merge(dict_df[['endCpG', 'end']], on=['endCpG'])
    df['end'] = df['end'] + 1
    df['endCpG'] = df['endCpG'] + 1
    df['chr'] = chrom
    return df[COORDS_COLS5]

def subset_by_chrom(cpgs, cf, chrom):
    chrom_ind = cf.chr.tolist().index(chrom)
    first = cf['size'][chrom_ind - 1] + 1 if chrom_ind > 0 else 1
    last = cf['size'][chrom_ind] + 1
    df_chrom = cpgs[(cpgs.startCpG>=first) & (cpgs.startCpG < last)]
    if (df_chrom.endCpG > last).any():
        eprint(f'[wt convert] ERROR: CpG range crosses chromosomes')
        eprint(df_chrom[df_chrom.endCpG > last].head(1).values.flatten())
        raise IllegalArgumentError('Invalid site file')
    return df_chrom

def add_bed_to_cpgs(cpgs, genome, threads):
    # load chromosomes sizes (in GpGs):
    cf = GenomeRefPaths(genome).get_chrom_cpg_size_table()
    cf['size'] = np.cumsum(cf['size'])

    params = []
    for chrom in sorted(set(cf.chr)):
        df_chrom = subset_by_chrom(cpgs, cf, chrom)
        if not df_chrom.empty:
            params.append((df_chrom.copy(), chrom, genome))
    assert len(params)
    p = Pool(threads)
    # arr = [chr_thread_cpg(*p) for p in params]
    arr = p.starmap(chr_thread_cpg, params)
    p.close()
    p.join()

    # concat chromosomes
    r = pd.concat(arr)
    # merge with original table, to keep the order
    r = cpgs.merge(r, how='left')[COORDS_COLS5]
    assert r.shape[0] == cpgs.shape[0]
    return r


####################################################################
#                                                                  #
#                      Add Annotations                             #
#                                                                  #
####################################################################

def get_anno(df, bed_path, genome):
    anno_path = GenomeRefPaths(genome).annotations
    if anno_path is None:
        return
    isnice, msg = is_block_file_nice(df)
    if not isnice:
        eprint(f'[wt convert] WARNING: No annotations added.\n'  \
               f'             Incompatible input BED file ({bed_path}):\n' \
               f'             {msg}')
        return
    if bed_path == sys.stdin:
        eprint(f'[wt convert] WARNING: No annotations added'  \
                f' when input file is passed with stdin.')
        return
    read = 'gunzip -cd' if bed_path.endswith('.gz') else 'cat'
    cmd = f'{read} {bed_path} | cut -f1-3 | sort -k1,1 -k2,2n | '
    cmd += f'bedtools intersect -wao -a - -b {anno_path} '
    cmd += '| bedtools merge -c 7,8 -o distinct,distinct '
    txt = subprocess.check_output(cmd, shell=True).decode()
    names = COORDS_COLS3 + ['type', 'gene']
    annodf = pd.read_csv(StringIO(txt), sep='\t', header=None, names=names)
    return annodf


def is_block_file_nice(df):

    # no duplicated blocks
    if (df.shape[0] != df.drop_duplicates().shape[0]):
        msg = 'Some blocks are duplicated'
        return False, msg

    # no overlaps between blocks
    sdf = df.sort_values(by='startCpG')
    if not (sdf['startCpG'][1:].values - sdf['endCpG'][:sdf.shape[0] - 1].values  >= 0).all():
        msg = 'Some blocks overlap'
        return False, msg

    return True, ''


####################################################################
#                                                                  #
#                            Main                                  #
#                                                                  #
####################################################################


def convert_single_region(args):
    gr = GenomicRegion(args)
    if args.parsable:
        r = gr.region_str if args.sites else '{}-{}'.format(*gr.sites)
    else:
        r = gr
    print(r)


def main():
    """
    Convert genomic region to CpG index range and vise versa
    """
    args = parse_args()

    try:
        if args.bed_file:
            convert_bed_file(args)
        elif args.site_file:
            convert_site_file(args)
        else:
            convert_single_region(args)

    except pd.errors.ParserError as e:
        print(f'Invalid input file.\n{e}')
        return


if __name__ == '__main__':
    main()

