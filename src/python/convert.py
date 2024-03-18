#!/usr/bin/python3 -u

import argparse
import os
import tempfile
from multiprocessing import Pool
import sys
import subprocess
import pandas as pd
import numpy as np
from utils_wgbs import add_GR_args, eprint, GenomeRefPaths, delete_or_skip, \
                       load_dict_section, add_multi_thread_args, IllegalArgumentError, \
                       read_shell, check_executable, COORDS_COLS3, COORDS_COLS5, \
                       add_loci_tool, validate_local_exe
from genomic_region import GenomicRegion


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    region_or_sites = add_GR_args(parser, bed_file=True, no_anno=True)
    region_or_sites.add_argument('--site_file',
                        help='text file with a single CpG indexes column,' \
                             ' or <startCpG, endCpG> columns.\n' \
                             'if "-" is passed, the file is read from stdin.')
    parser.add_argument('--out_path', '-o', help='Output path for bed file [stdout]')
    parser.add_argument('-d', '--debug', action='store_true')
    # parser.add_argument('--bedtools', action='store_true')
    parser.add_argument('-p', '--parsable', action='store_true',
                        help='Output a parsing friendly format')
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
    out_path = sys.stdout if args.out_path is None else args.out_path
    if not delete_or_skip(out_path, args.force):
        return
    # add CpG columns
    bed_file = args.bed_file
    # TODO: support stdin for -L in all wgbstools features, and add it to the help message
    if bed_file == '-':
        bed_file = sys.stdin
    add_anno = (not args.parsable) and (not args.no_anno)

    if not check_executable('bedtools', verbose=False):
        # eprint('continue with a slower implementation')
        r = add_cpgs_to_bed(bed_file=bed_file,
                genome=args.genome,
                drop_empty=args.drop_empty,
                threads=args.threads,
                add_anno=add_anno)
    else:
        r = bedtools_conversion(bed_file, args.genome,
                args.drop_empty, add_anno, args.debug)
    r.to_csv(out_path, sep='\t', header=None, index=None, na_rep='NA')


def load_bed(bed_path, nrows=None):
    try:
        # TODO: handling a bed with a header line. Test
        df = pd.read_csv(bed_path, sep='\t', header=None, nrows=nrows, comment='#')
        df.columns = COORDS_COLS3 + list(df.columns)[3:]
        if not (str(df.iloc[0, 1]).isdigit() and str(df.iloc[0, 2]).isdigit()):
            eprint('[wt convert] Header line detected. Ignoring first line of input')
            df = df.iloc[1: , :].reset_index(drop=True)
            df = df.astype({'start': int, 'end': int})
        return df
    except pd.errors.EmptyDataError:
        eprint('[wt convert] ERROR: empty bed file')
        raise IllegalArgumentError('Invalid bed file')

# TODO: in some cases it differs from the slow implementations
#       Try -L /cs/cbio/tommy/indexes/hg19.CDS.bed
# TODO2: It fails for >3 columns (bedtools itself fails)
# TODO3: implement my own c++ convert tool
def bedtools_conversion(bed_file, genome, drop_empty, add_anno, debug):
    df = load_bed(bed_file)
    tmp_name = tempfile.NamedTemporaryFile().name
    df.sort_values(by=['chr', 'start', 'end']).iloc[:, :3].\
        drop_duplicates().to_csv(tmp_name, sep='\t', header=None, index=None)
    ref = GenomeRefPaths(genome).dict_path
    cmd = f"tabix -R {tmp_name} {ref} | "
    cmd += "awk -v OFS='\t' '{print $1,$2,$2+1,$3}' | "
    cmd += " sort -k1,1 -k2,2n -u | "       # sort is required for cases of overlapping blocks
    cmd += f"bedtools intersect -sorted -b - -a {tmp_name} -loj | "
    cmd += "bedtools groupby -g 1,2,3 -c 7,7 -o first,last | "
    cmd += " awk -v OFS='\t' '{print $1,$2,$3,$4,$5+1;}' "
    cmd += "| sed 's/\.\t1/NA\tNA/g'"       # replace missing values with (NA)s
    if debug:
        eprint(cmd.replace('\t', '\\t'))
    rf = read_shell(cmd, names=COORDS_COLS5)
    if not debug:
        os.unlink(tmp_name)

    # if there are missing values, the CpG columns' type
    # will be float or object. Change it to Int64
    if rf.empty:
        raise IllegalArgumentError('[wt convert] Error: failed with bedtools wrapping')

    if rf['startCpG'].dtype != int:
        rf = rf.astype({'startCpG': 'Int64', 'endCpG': 'Int64'})

    df = df.merge(rf, how='left', on=COORDS_COLS3)
    df = df[COORDS_COLS5 + list(df.columns)[3:-2]]

    if drop_empty:
        df.dropna(inplace=True, subset=['startCpG', 'endCpG'])

    # add annotations:
    if add_anno:
        df = get_anno(df, genome, bed_file)

    return df

def slow_conversion(df, genome):
    df = df.iloc[:, :3]
    startCpGs = []
    endCpGs = []
    for _, row in df.iterrows():
        try:
            sites = GenomicRegion(region='{}:{}-{}'.format(*row), genome_name=genome).sites
        except IllegalArgumentError:
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
                    'Found overlaps in the input bed file. Conversion may be slow.\n' \
                    '             Install bedtools for better performance')
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
        r = get_anno(r, genome, bed_file)

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
    out_path = sys.stdout if args.out_path is None else args.out_path
    if not delete_or_skip(out_path, args.force):
        return
    # add loci columns
    add_bed_to_cpgs(args.site_file, args.genome, args.out_path)


def add_bed_to_cpgs(site_file, genome, out_path=None):
    validate_local_exe(add_loci_tool)
    g = GenomeRefPaths(genome)
    cmd = f'cat {site_file} | {add_loci_tool} {g.dict_path} {g.chrom_cpg_sizes}'
    if (out_path is not None) and out_path != sys.stdout:
        cmd += f' > {out_path}'
    subprocess.check_call(cmd, shell=True)


####################################################################
#                                                                  #
#                      Add Annotations                             #
#                                                                  #
####################################################################

def get_anno(df, genome, bed_file):
    try:
        anno_path = GenomeRefPaths(genome).annotations
        if anno_path is None:
            return df
        cmd = 'gunzip -c' if bed_file.endswith('gz') else 'cat'
        cmd += f' {bed_file} | cut -f1-3 | sort -k1,1 -k2,2n -u | '
        cmd += f'bedtools intersect -a - -b {anno_path} -wao | '
        cmd += 'bedtools merge -i - -c 7,8 -o distinct,distinct'
        names = COORDS_COLS3 + ['type', 'gene']
        rf = read_shell(cmd, names=names)
        return df.merge(rf, how='left', on=COORDS_COLS3)
    except Exception as e:
        eprint('[wt convert] WARNING: No annotations added.')
        eprint(cmd)
        eprint(e)
        return df


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
