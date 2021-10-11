#!/usr/bin/python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from utils_wgbs import load_beta_data, validate_file_list, add_GR_args, eprint
from genomic_region import GenomicRegion
import os.path as op
from itertools import combinations


def comp2(a, b, ids, cov_thresh, out_dir):
    # Ignore sites with coverage lower than cov_thresh in at
    # least one of the two compared files
    high_cov = np.min(np.c_[a[:, 1], b[:, 1]], axis=1) >= cov_thresh
    def table2vec(table):
        table = table[high_cov]
        return table[:, 0] / table[:, 1]
    plt.hist2d(table2vec(a), table2vec(b), bins=20, cmap=plt.cm.jet, norm=LogNorm())
    plt.title('{} : {}'.format(*ids), fontsize=10)
    if out_dir is not None:
        outpath = op.join(out_dir, '{}:{}.png'.format(*ids))
        plt.savefig(outpath)
        eprint(f'[wt cmp] dumped figure to {outpath}')


# def compare_all_paires(betas, min_cov, sites):
def compare_all_paires(args):
    betas = args.betas
    sites = GenomicRegion(args).sites
    tables = [load_beta_data(b, sites) for b in betas]
    names = [op.splitext(op.basename(b))[0] for b in betas]

    for x, y in combinations(range(len(tables)), r=2):
        plt.figure()
        comp2(tables[x], tables[y], (names[x], names[y]), args.min_cov, args.out_dir)
    if args.show or args.out_dir is None:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('betas', nargs='+')
    parser.add_argument('--out_dir', '-o',
            help='Dump figures to this directory. If not specified, --show flag is set')
    parser.add_argument('--show', action='store_true',
            help='Display the figures using matplotlib.pyplot.show.')
    parser.add_argument('--min_cov', '-c', type=int, default=10,
                        help='Minimal coverage to consider. '
                             'Sites with coverage lower than this value are ignored')
    add_GR_args(parser)
    args = parser.parse_args()
    return args


def main():
    """
    Compare between pairs of beta files, by plotting a 2d histogram
    for every pair.
    Drop sites with low coverage (< cov_thresh argument),
    for performance and robustness.
    """
    args = parse_args()
    validate_file_list(args.betas, '.beta', min_len=2)
    # compare_all_paires(args.betas, args.min_cov, GenomicRegion(args).sites)
    compare_all_paires(args)


if __name__ == '__main__':
    main()
