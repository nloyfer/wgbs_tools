#!/usr/bin/python3

import numpy as np
import os.path as op
import sys
import os
import pandas as pd
import argparse
from utils_wgbs import load_beta_data, add_GR_args, \
        add_multi_thread_args, load_dict_section
from multiprocessing import Pool
from genomic_region import GenomicRegion
from beta_cov import pretty_name
import tempfile
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('betas', nargs='+', help='one or more beta files')
    parser.add_argument('--width', '-w', type=int, default=120, help='max width to print output table [120]')
    add_GR_args(parser, bed_file=True)
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args


def print_stats(beta_path, data):
    assert data.size, beta_path + ': Data table is empty!'
    name = pretty_name(beta_path)#[:20]
    names = []
    vals = []
    np.seterr(divide='ignore', invalid='ignore')

    vals.append(str(np.nanmean(data[:, 0] / data[:, 1] * 100).round(2)))
    names.append('mean meth. (%)')

    vals.append(f'{(data[:, 1] > 0).sum():,}')
    names.append('covered sites')

    # vals.append(f'{(data[:, 1] >= 5).sum():,}')
    # names.append('covered sites (5+)')

    vals.append(f'{(data[:, 1] >= 10).sum():,}')
    names.append('covered sites (10+)')

    vals.append(f'{(data[:, 1]).max():,}')
    names.append('max depth')

    vals.append(f'{(data[:, 1]).mean().round(2):,}')
    names.append('mean depth')

    # print results
    df = pd.DataFrame(data={'names': names, name: vals})
    df = df.set_index('names')
    return df


def load_beta_by_bed(beta_path, bed_path, genome_name):
    tmp = tempfile.NamedTemporaryFile(delete=False)
    data = np.array(())
    try:
        cmd = f'cut -f 1-3 {bed_path} > {tmp.name}'
        subprocess.check_call(cmd, shell=True)
        inds = load_dict_section(' -R ' + tmp.name, genome_name)['idx'].values - 1
        data = load_beta_data(beta_path)[np.unique(inds), :]
    finally:
        tmp.close()
        os.unlink(tmp.name)
    return data


def beta_stats(beta_path, gr=None, bed_path=None):
    if bed_path:
        data = load_beta_by_bed(beta_path, bed_path, gr.genome.genome)
    else:
        data = load_beta_data(beta_path, gr.sites, gr.genome)
    return print_stats(beta_path, data)


def main():
    """
    Print global stats of one or more beta/lbeta file(s)
    """

    args = parse_args()

    gr = GenomicRegion(args)

    params = [(beta, gr, args.bed_file) for beta in args.betas]
    p = Pool(args.threads)
    stats_arr = p.starmap(beta_stats, params)
    p.close()
    p.join()

    df = pd.concat(stats_arr, axis=1)
    pd.set_option('display.max_columns', None)  # Display all columns
    pd.set_option('display.max_rows', None)     # Display all rows
    pd.set_option('display.width', args.width)  # Set max width
    print(df.T)


if __name__ == '__main__':
    main()

