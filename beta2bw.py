#!/usr/bin/python3 -u

import argparse
import os.path as op
import sys
import subprocess
import numpy as np
from utils_wgbs import delete_or_skip, load_beta_data, validate_files_list, load_dict, GenomeRefPaths
import os


BG_EXT = '.bedGraph'
BW_EXT = '.bw'
COV_BG_EXT = '_cov' + BG_EXT
COV_BW_EXT = '_cov' + BW_EXT


def bedGraph2bigwig(bedGraph_path, keepBG, chrom_sizes):

    # Sort the bedGraph
    cmd = 'sort -k1,1 -k2,2n -o {b} {b}'.format(b=bedGraph_path)
    r = subprocess.call(cmd, shell=True)
    if r:
        print('Failed sorting bedGraph.\nFile: {}\nError code: {}'.format(bedGraph_path, r), file=sys.stderr)
        print('cmd', cmd, file=sys.stderr)
        return

    # Convert bedGraph to bigWig:
    bw_path = bedGraph_path.replace(BG_EXT, BW_EXT)
    r = subprocess.call(['bedGraphToBigWig', bedGraph_path, chrom_sizes, bw_path])
    if r:
        print('Failed converting bedGraph to bigWig. File: {}\nError code: {}'.format(bedGraph_path, r),
              file=sys.stderr)
        return

    # pigz or delete the bedGraph:
    if keepBG:
        subprocess.call(['pigz', '-f', bedGraph_path])
    else:
        os.remove(bedGraph_path)


def betas2bws(args):
    chrom_sizes = GenomeRefPaths(args.genome).chrom_sizes

    rf = None       # Reference dictionary
    for beta in args.beta_paths:
        go = True
        print('Converting {}...'.format(op.basename(beta)), file=sys.stderr)
        # Check if beta should be skipped:
        outpath = op.join(args.outdir, op.splitext(op.basename(beta))[0])
        for suf in (BW_EXT, COV_BW_EXT):
            if not delete_or_skip(outpath + suf, args.force):
                go = False

        if not go:
            continue

        # Load dict (at most once) and beta
        if rf is None:
            rf = load_dict(genome_name=args.genome)
            rf['end'] = rf['start'] + 1
        barr = load_beta_data(beta)

        # dump coverage:
        rf['cov'] = barr[:, 1]
        rf[rf['cov'] > 0].to_csv(outpath + COV_BG_EXT, sep='\t', header=None, index=None)
        del rf['cov']

        # dump beta values
        # betav = barr[:, 0] / barr[:, 1].astype(np.float16)
        betav = np.divide(barr[:, 0], barr[:, 1].astype(np.float16), where=barr[:, 1] > 0)
        betav[barr[:, 1] == 0] = -1
        rf['beta'] = betav
        rf.to_csv(outpath + BG_EXT, sep='\t', header=None, index=None)
        del rf['beta']

        # convert bedGraphs to bigWigs:
        bedGraph2bigwig(outpath + COV_BG_EXT, args.bedGraph, chrom_sizes)
        bedGraph2bigwig(outpath + BG_EXT, args.bedGraph, chrom_sizes)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('beta_paths', nargs='+')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-b', '--bedGraph', action='store_true', help='Keep bedGraphs as well as bigwigs')
    parser.add_argument('--outdir', '-o', default='.')
    parser.add_argument('--genome', help='Genome reference name. Default is hg19.', default='hg19')
    args = parser.parse_args()
    return args


def main():
    """
    Convert beta file[s] to bigbiw file[s].
    """
    args = parse_args()
    validate_files_list(args.beta_paths, '.beta')
    betas2bws(args)


if __name__ == '__main__':
    main()

