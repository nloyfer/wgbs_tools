#!/usr/bin/python3

import os.path as op
import argparse
from multiprocessing import Pool
import numpy as np
from utils_wgbs import load_beta_data, add_GR_args, \
        add_multi_thread_args
from genomic_region import GenomicRegion
from beta_to_blocks import collapse_process, is_block_file_nice, load_blocks_file


def plot_hist2(beta_path, sites):
    ymax = 255 if beta_path.endswith('.beta') else 1000
    import plotille
    data = load_beta_data(beta_path, sites)[:, 1]
    h = plotille.histogram(data, x_min=0, x_max=ymax)
    print(h)

def plot_hist(names, covs):
    import matplotlib.pyplot as plt
    plt.rcdefaults()
    plt.hist(covs)
    plt.title(f'beta coverage histogram\nmean cov:{np.mean(covs):.2f}')

    plt.figure()
    y_pos = np.arange(len(covs))
    plt.bar(y_pos, covs)
    plt.xticks(y_pos, names, rotation='vertical')
    plt.subplots_adjust(bottom=0.15)
    plt.ylabel('Coverage')
    plt.title('Coverage bar chart')
    plt.tight_layout()
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('betas', nargs='+', help='one or more beta files')
    parser.add_argument('--plot', action='store_true', help='Plot histogram of coverages')
    parser.add_argument('--hist', action='store_true', help='Plot in-terminal histogram of coverages')
    add_GR_args(parser, bed_file=True)
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args


def pretty_name(beta_path):
    return op.splitext(op.basename(beta_path))[0]


def beta_cov_by_bed(beta_path, blocks_df):
    nr_sites_covered = (blocks_df['endCpG'] - blocks_df['startCpG']).sum()
    if not nr_sites_covered:
        return 0
    is_nice = is_block_file_nice(blocks_df)[0]
    data = collapse_process(beta_path, blocks_df, is_nice)
    total_cov = data[:, 1].sum()
    return total_cov / nr_sites_covered


def beta_cov(beta_path, sites=None, blocks_df=None, print_res=False):
    if blocks_df is not None:
        res = beta_cov_by_bed(beta_path, blocks_df)
    else:
        res = np.mean(load_beta_data(beta_path, sites)[:, 1])
    if print_res:
        print('{}\t{:.2f}'.format(pretty_name(beta_path), res))
    return res


def main():
    """
    Calculate the average coverage of one or more beta files.
    Print the results.
    """

    args = parse_args()


    sites = GenomicRegion(args).sites

    blocks_df = load_blocks_file(args.bed_file) if args.bed_file else None

    params = [(beta, sites, blocks_df, False) for beta in args.betas]
    # covs = [beta_cov(*p) for p in params]
    # return
    p = Pool(args.threads)
    covs = p.starmap(beta_cov, params)
    p.close()
    p.join()

    for cov, beta_path in zip(covs, args.betas):
        print('{}\t{:.2f}'.format(pretty_name(beta_path), cov))
        if args.hist:
            plot_hist2(beta_path, sites)
    if args.plot:
        plot_hist([pretty_name(b) for b in args.betas], covs)


if __name__ == '__main__':
    main()
