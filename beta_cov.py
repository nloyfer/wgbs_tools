#!/usr/bin/python3

import numpy as np
import os.path as op
import argparse
from utils_wgbs import load_beta_data, add_GR_args, BedFileWrap
import multiprocessing
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


def pretty_name(beta_path):
    return op.splitext(op.basename(beta_path))[0]


def beta_cov_by_bed(beta_path, bed_wrapper):
    nr_sites = 0
    total_cov = 0
    for gr in bed_wrapper.iter_grs():
        table = load_beta_data(beta_path, gr.sites)[:, 1]
        nr_sites += table.size
        total_cov += table.sum()
    return total_cov / nr_sites if nr_sites else 0


def beta_cov(beta_path, sites=None, bed_wrapper=None, print_res=False):
    if bed_wrapper:
        res = beta_cov_by_bed(beta_path, bed_wrapper)
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

    bedw = BedFileWrap(args.bed_file) if args.bed_file else None

    processes = []
    with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
        for beta in args.betas:
            params = (beta, GenomicRegion(args).sites, bedw, True)
            processes.append(p.apply_async(beta_cov, params))
        p.close()
        p.join()
    covs = [pr.get() for pr in processes]

    if args.plot:
        plot_hist([pretty_name(b) for b in args.betas], covs)


if __name__ == '__main__':
    main()
