#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_files_list, add_GR_args, MAX_READ_LEN, splitextgz, BedFileWrap, eprint
from genomic_region import GenomicRegion
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool
import sys
import os.path as op


def awk_wrap(awk_cmd):
    x = subprocess.check_output(awk_cmd, shell=True).decode()
    return np.array([int(r) if r else 0 for r in x.split('\n')[:-1]], dtype=int, ndmin=2)


class FragLen:
    def __init__(self, unq, args, gr=None):
        self.args = args
        self.unq = unq
        self.gr = gr if gr else GenomicRegion(args)
        m = args.max_frag_size
        self.fill_arr_cmd = ' {if ($3 > %s) {$3 = %s}; arr[$3] += $5}' % (m, m)
        self.print_arr_cmd = ' END {for (x=1; x <= %s; x++) print arr[x]}\'' % m

    def run_small_region(self):
        start, end = self.gr.bp_tuple
        cmd = 'tabix {} '.format(self.unq)
        cmd += '{}:{}-{} '.format(self.gr.chrom, max(1, int(start) - MAX_READ_LEN), end)
        cmd += ' | awk \'{if (($2 + $3) > %s) ' % start
        cmd += self.fill_arr_cmd + '}' + self.print_arr_cmd
        return awk_wrap(cmd)

    def chrom_cmd(self, chrom):
        cmd = 'tabix {} {} | awk \''.format(self.unq, chrom)
        cmd = cmd + self.fill_arr_cmd + self.print_arr_cmd
        return cmd

    def run_whole_genome(self):
        chroms = list(pd.read_csv(self.gr.genome.chrom_sizes, sep='\t', header=None, usecols=[0])[0])
        with Pool() as p:
            pr = [p.apply_async(awk_wrap, (self.chrom_cmd(c),)) for c in chroms]
            p.close()
            p.join()
        return np.sum(np.concatenate([r.get() for r in pr]), axis=0)


def compose_fig_path(unq, outdir, grs):
    if not outdir:
        return
    res = op.join(outdir, op.basename(splitextgz(unq)[0]))
    if grs and len(grs) == 1:
        res += '.{}'.format(grs[0].region_str)
    res += '.png'
    return res


def run_single_unq(unq, grs, args):
    eprint(unq)
    if not grs:  # process the whole unq file (no -L,-s,-r was specified)
        x = FragLen(unq, args).run_whole_genome()
    else:
        x = np.sum([FragLen(unq, args, gr).run_small_region() for gr in grs], axis=0)

    # print values to stdout:
    if args.verbose:
        np.savetxt(sys.stdout, x.reshape((1, -1)), fmt='%s', delimiter=' ')

    # plot:
    plt.figure()
    plt.plot(np.arange(x.size), x.flatten())
    plt.title('Fragment lengths\n' + op.basename(splitextgz(unq)[0]))

    # dump figure:
    if args.outdir:
        plt.savefig(compose_fig_path(unq, args.outdir, grs))


def multi_FragLen(args):
    if args.bed_file and (args.region or args.sites):
        eprint('-L, -s and -r are mutually exclusive')
        return

    if args.region or args.sites:
        grs = [GenomicRegion(args)]
    elif args.bed_file:
        grs = BedFileWrap(args.bed_file).iter_grs()
    else:
        grs = []

    for unq in args.unq_paths:
        run_single_unq(unq, grs, args)

    if args.display:
        plt.show()


##########################
#                        #
#         Main           #
#                        #
##########################


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('unq_paths', nargs='+')
    parser.add_argument('-v', '--verbose', action='store_true', help='print the histogram values to stdout')
    parser.add_argument('-m', '--max_frag_size', type=int, default=500,
                        help='Maximum fragment size. Longer fragments will be trimmed. [500]')
    parser.add_argument('--outdir', '-o', help='output directory for the histogram figure(s) [None]')
    parser.add_argument('--display', action='store_true', help='Display histogram plot(s) (plt.show)')
    parser.add_argument('-L', '--bed_file',
                        help='pat: Only output reads overlapping the input BED FILE')
    add_GR_args(parser)
    return parser.parse_args()


def main():
    """
    Plot histogram of reads lengths of unq file
    Output to stdout the histogram values if requested
    """
    args = parse_args()
    validate_files_list(args.unq_paths, 'unq.gz')
    multi_FragLen(args)


if __name__ == '__main__':
    main()
