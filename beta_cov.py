#!/usr/bin/python3

import numpy as np
import os.path as op
import argparse
from utils_wgbs import load_beta_data, add_GR_args, BedFileWrap
from genomic_region import GenomicRegion


def plot_hist(names, covs):
    import matplotlib.pyplot as plt
    plt.rcdefaults()
    plt.hist(covs)
    plt.title('beta coverage histogram\nmean cov:{:.2f}'.format(np.mean(covs)))

    plt.figure()
    y_pos = np.arange(len(covs))
    plt.bar(y_pos, covs)
    plt.xticks(y_pos, names, rotation='vertical')
    plt.subplots_adjust(bottom=0.15)
    plt.ylabel('Coverage')
    plt.title('Coverage bar chart')
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('betas', nargs='+', help='one or more beta files')
    parser.add_argument('--plot', action='store_true', help='Plot histogram of coverages')
    add_GR_args(parser)
    parser.add_argument('-L', '--bed_file', help='Only output coverage overlapping the input BED FILE. ')
    args = parser.parse_args()
    return args


def beta_cov_by_bed(beta_path, bed_wrapper):
    nr_sites = 0
    total_cov = 0
    for gr in bed_wrapper.iter_grs():
        table = load_beta_data(beta_path, gr.sites)[:, 1]
        nr_sites += table.size
        total_cov += table.sum()
    return total_cov / nr_sites


def beta_cov(beta_path, sites=None, bed_wrapper=None):
    if bed_wrapper:
        return beta_cov_by_bed(beta_path, bed_wrapper)
    else:
        return np.mean(load_beta_data(beta_path, sites)[:, 1])


def main():
    """
    Calculate the average coverage of one or more beta files.
    Print the results.
    """
    args = parse_args()

    bedw = BedFileWrap(args.bed_file) if args.bed_file else None
    covs = []
    names = []

    for beta in args.betas:
        cov = beta_cov(beta, sites=GenomicRegion(args).sites, bed_wrapper=bedw)
        covs.append(cov)
        name = op.splitext(op.basename(beta))[0]
        names.append(name)
        print('{}\t{:.2f}'.format(name, cov))

    if args.plot:
        plot_hist(names, covs)


if __name__ == '__main__':
    main()

