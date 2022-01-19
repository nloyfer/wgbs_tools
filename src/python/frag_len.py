#!/usr/bin/python3 -u

import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool
import sys
import os.path as op
from utils_wgbs import validate_file_list, add_GR_args, splitextgz, eprint, main_script
from genomic_region import GenomicRegion


def awk_wrap(awk_cmd):
    # print(awk_cmd)
    x = subprocess.check_output(awk_cmd, shell=True).decode()
    return np.array([int(r) if r else 0 for r in x.split('\n')[:-1]], dtype=int, ndmin=2)


class FragLen:
    def __init__(self, pat, args, gr=None):
        self.args = args
        self.pat = pat
        self.gr = gr if gr else GenomicRegion(args)
        m = args.max_frag_size
        self.awk_hist_cmd = ' | awk \'{l=length($3) ; if (l > %s) {l = %s}; arr[l] += $4} END {for (x=1; x <= %s; x++) print arr[x]}\'' % (m, m, m)

    def run_small_region(self):
        cmd = f'{main_script} cview {self.pat} -r {self.gr.region_str} {self.awk_hist_cmd}'
        return awk_wrap(cmd)

    def run_bed(self):
        cmd = f'{main_script} cview {self.pat} -L {self.args.bed_file} {self.awk_hist_cmd}'
        return awk_wrap(cmd)

    def chrom_cmd(self, chrom):
        return f'tabix {self.pat} {chrom} {self.awk_hist_cmd}'

    def run_whole_genome(self):
        chroms = self.gr.genome.get_chroms()
        with Pool() as p:
            pr = [p.apply_async(awk_wrap, (self.chrom_cmd(c),)) for c in chroms]
            p.close()
            p.join()
        return np.sum(np.concatenate([r.get() for r in pr]), axis=0)


def compose_fig_path(pat, outdir):
    if outdir:
        return op.join(outdir, op.basename(splitextgz(pat)[0])) + '.png'


def plot_hist(data, max_frag_size, pat):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.arange(1, data.size + 1), data)
    major_ticks = np.arange(1, max_frag_size, 5)
    minor_ticks = np.arange(1, max_frag_size, 1)

    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)

    # Or if you want different settings for the grids:
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)
    plt.ylim(bottom=0)
    plt.xlim(left=1)
    plt.title('Fragment lengths (CpGs)\n' + op.basename(splitextgz(pat)[0]))


def run_single_pat(pat, args):
    eprint(pat)
    fl = FragLen(pat, args)
    if args.region or args.sites:
        x = fl.run_small_region()
    elif args.bed_file:
        x = fl.run_bed()
    else:
        x = fl.run_whole_genome()

    if not x.sum():
        eprint(f'[wt frag] Empty list of lengths for {pat}')
        return

    # print values to stdout:
    if args.verbose:
        np.savetxt(sys.stdout, x.reshape((1, -1)), fmt='%s', delimiter=' ')

    # plot:
    if args.outdir or args.display:
        if args.verbose:
            eprint('[wt frag] plotting...')
        plot_hist(x.flatten(), args.max_frag_size, pat)

    # dump figure:
    if args.outdir:
        fpath = compose_fig_path(pat, args.outdir)
        if args.verbose:
            eprint(f'[wt frag] dumping {fpath}...')
        plt.savefig(fpath)


def multi_FragLen(args):
    for pat in args.pat_paths:
        run_single_pat(pat, args)

    if args.display:
        plt.show()


##########################
#                        #
#         Main           #
#                        #
##########################


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat_paths', nargs='+')
    parser.add_argument('-v', '--verbose', action='store_true', help='print the histogram values to stdout')
    parser.add_argument('-m', '--max_frag_size', type=int, default=30,
                        help='Maximum fragment size. Longer fragments will be trimmed. [30]')
    parser.add_argument('--outdir', '-o', help='output directory for the histogram figure(s) [None]')
    parser.add_argument('--display', action='store_true', help='Display histogram plot(s) (plt.show)')
    add_GR_args(parser, bed_file=True)
    return parser.parse_args()


def main():
    """
    Plot histogram of reads lengths (in sites) of pat file
    Output to stdout the histogram values if requested
    """
    args = parse_args()
    validate_file_list(args.pat_paths, '.pat.gz')
    multi_FragLen(args)


if __name__ == '__main__':
    main()
