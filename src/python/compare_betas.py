#!/usr/bin/python3
import argparse
import os.path as op
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from utils_wgbs import load_beta_data, validate_file_list, add_GR_args, eprint
from genomic_region import GenomicRegion


def comp2(a, b, cov_thresh, bins, ax):
    # Ignore sites with coverage lower than cov_thresh in at
    # least one of the two compared files
    high_cov = np.min(np.c_[a[:, 1], b[:, 1]], axis=1) >= cov_thresh
    def table2vec(table):
        table = table[high_cov]
        return table[:, 0] / table[:, 1]
    ax.hist2d(table2vec(b), table2vec(a), bins=bins, cmap=plt.cm.jet, norm=LogNorm())
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 1)


def compare_all_paires(args):
    # breakpoint()
    betas = args.betas
    sites = GenomicRegion(args).sites
    tables = [load_beta_data(b, sites) for b in betas]
    names = [op.splitext(op.basename(b))[0] for b in betas]
    # break names to lines
    nnames = []
    k = 20
    for n in names:
        lst = [n[0+i:k+i] for i in range(0, len(n), k)]
        nn = '\n'.join(lst)
        nnames.append(nn)

    N = len(tables)
    fig, axs = plt.subplots(N, N)
    for i in range(N):
        for j in range(i + 1):
            comp2(tables[i], tables[j], args.min_cov, args.bins, axs[i, j])
        axs[i, 0].set_ylabel(nnames[i], fontsize=8)
    for j in range(N):
        axs[N - 1, j].set_xlabel(nnames[j], fontsize=8)

    # remove empty subplots
    for i in range(N):
        for j in range(i + 1, N):
            fig.delaxes(axs[i, j])

    for ax in axs.flat:
        ax.label_outer()

    fig.tight_layout()


    if args.outpath is not None:
        plt.savefig(args.outpath)
        eprint(f'[wt cmp] dumped figure to {args.outpath}')

    if args.show or args.outpath is None:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('betas', nargs='+')
    parser.add_argument('--outpath', '-o',
            help='Dump figure to this path (e.g., pdf/png). If not specified, --show flag is set')
    parser.add_argument('--show', action='store_true',
            help='Display the figures using matplotlib.pyplot.show.')
    parser.add_argument('--min_cov', '-c', type=int, default=10,
                        help='Minimal coverage to consider. '
                             'Sites with coverage lower than this value are ignored')
    parser.add_argument('--bins', type=int, default=101,
                        help='Histogram bins (resolution) [101]')
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
    validate_file_list(args.betas, min_len=2)
    compare_all_paires(args)


if __name__ == '__main__':
    main()
