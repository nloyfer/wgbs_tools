#!/usr/bin/python3

import numpy as np
import os.path as op
import sys
import pandas as pd
import argparse
from utils_wgbs import load_beta_data, add_GR_args, \
        add_multi_thread_args
from multiprocessing import Pool
from genomic_region import GenomicRegion
from beta_cov import pretty_name


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('betas', nargs='+', help='one or more beta files')
    parser.add_argument('--width', '-w', type=int, default=120, help='max width to print output table [120]')
    add_GR_args(parser, bed_file=False)  # TODO: support bed
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args


def print_stats(fpath, sites):
    data = load_beta_data(fpath, sites)
    name = pretty_name(fpath)#[:20]
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


def beta_stats(beta_path, sites=None):
    return print_stats(beta_path, sites)


def main():
    """
    Print global stats of one or more beta/lbeta file(s)
    """

    args = parse_args()

    sites = GenomicRegion(args).sites

    params = [(beta, sites) for beta in args.betas]
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

