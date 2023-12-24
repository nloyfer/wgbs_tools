#!/usr/bin/python3 -u

import argparse
import os.path as op
import sys
import subprocess
from utils_wgbs import delete_or_skip, validate_file_list, GenomeRefPaths, \
                       add_GR_args, IllegalArgumentError, check_executable, \
                       safe_remove
from genomic_region import GenomicRegion
from beta2bed import beta_to_bed, beta2bed_build_cmd
from cview import subprocess_wrap_sigpipe
import os
import numpy as np
from pathlib import Path

BG_EXT = '.bedGraph'
BW_EXT = '.bigwig'

def b2bw_log(*args, **kwargs):
    print('[ wt beta_2bw ]', *args, file=sys.stderr, **kwargs)


class BetaToBigWig:
    def __init__(self, args):
        self.args = args
        self.gr = GenomicRegion(args)
        self.outdir = args.outdir
        self.name = ''
        if not op.isdir(self.outdir):
            raise IllegalArgumentError('Invalid output directory: ' + self.outdir)
        self.chrom_sizes = GenomeRefPaths(args.genome).chrom_sizes

    def bed_graph_to_bigwig(self, bed_graph, bigwig, is_sorted=True):
        """
        Generate a bigwig file from a bedGraph
        :param bed_graph: path to bedGraph (input)
        :param bigwig: path to bigwig (output)
        """

        # Convert bedGraph to bigWig:
        if not is_sorted:
            b2bw_log(f'[{self.name}] sort bedgraph...')
            subprocess.check_call('sort -k1,1 -k2,2n {o} -o {o}'.format(o=bed_graph), shell=True)
        b2bw_log(f'[{self.name}] convert bed to bigwig...')
        subprocess.check_call(['bedGraphToBigWig', bed_graph, self.chrom_sizes, bigwig])

        # compress or delete the bedGraph:
        if self.args.bedGraph:
            compress = 'pigz' if check_executable('pigz') else 'gzip'
            subprocess.check_call([compress, '-f', bed_graph])
        else:
            os.remove(bed_graph)


    def run_beta_to_bw(self, beta_path):
        self.name = op.basename(beta_path)

        prefix = op.join(self.outdir, op.splitext(self.name)[0])
        out_bigwig = prefix + BW_EXT
        out_bed_graph = prefix + BG_EXT

        # Check if the current file should be skipped:
        if not delete_or_skip(out_bigwig, self.args.force):
            return

        # convert beta to bed:
        if self.gr.is_whole():
            # sort chromosomes to unix/ucsc order (chr1, chr10, chr11, ..., chr2, ...)
            cmd = f'sort -k1,1 {self.chrom_sizes} | cut -f1'
            sorted_chroms = subprocess.check_output(cmd, shell=True).decode().split()

            # dump bedGraph chromosome by chromosome, so there is no need to sort
            safe_remove(out_bed_graph)
            Path(out_bed_graph).touch()
            b2bw_log(f'[{self.name}] Dumping bedGraph...')
            for chrom in sorted_chroms:
                # b2bw_log(f'[{self.name}] {chrom}')
                cmd = beta2bed_build_cmd(beta_path=beta_path,
                                   gr=GenomicRegion(region=chrom, genome_name=self.gr.genome_name),
                                   min_cov=self.args.min_cov,
                                   bed_file=None,
                                   mean=True,
                                   keep_na=self.args.keep_na)
                cmd += f' >> {out_bed_graph}'
                subprocess_wrap_sigpipe(cmd)
        else:
            beta_to_bed(beta_path=beta_path,
                        gr=self.gr,
                        bed_file=self.args.bed_file,
                        min_cov=self.args.min_cov,
                        mean=True,
                        keep_na=self.args.keep_na,
                        force=True,
                        opath=out_bed_graph)

        # convert bedGraphs to bigWigs:
        self.bed_graph_to_bigwig(out_bed_graph, out_bigwig, is_sorted=(not self.args.bed_file))

        if self.args.dump_cov:
            out_cov_bigwig = prefix + '.cov' + BW_EXT
            out_cov_bed_graph = prefix + '.cov' + BG_EXT
            # convert beta to bed:
            b2bw_log(f'[{self.name}] Dumping cov bed...')
            cmd = beta2bed_build_cmd(beta_path=beta_path,
                               gr=self.gr,
                               bed_file=self.args.bed_file,
                               min_cov=self.args.min_cov,
                               mean=False,
                               keep_na=self.args.keep_na)
            cmd += " | awk '{print $1,$2,$3,$5}' | sort -k1,1 -k2,2n"
            cmd += f' > {out_cov_bed_graph}'
            subprocess_wrap_sigpipe(cmd)

            # convert bedGraphs to bigWigs:
            self.bed_graph_to_bigwig(out_cov_bed_graph, out_cov_bigwig)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('beta_paths', nargs='+')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Overwrite existing files if existed')
    parser.add_argument('--keep_na', action='store_true',
                        help='If set, missing CpG sites are not removed from the output'
                             ' They are assigned with a "-1" value.')
    parser.add_argument('-b', '--bedGraph', action='store_true',
                        help='Keep (gzipped) bedGraphs as well as bigwigs')
    parser.add_argument('--dump_cov', action='store_true',
                        help='Generate coverage bigiwig in addition to beta values bigwig')
    parser.add_argument('-c', '--min_cov', type=int, default=1,
                        help='Minimal coverage to consider when computing beta values.'
                             ' Default is 1 (include all observations). '
                             ' Sites with less than MIN_COV coverage are considered as missing.')
    parser.add_argument('--outdir', '-o', default='.', help='Output directory. [.]')
    add_GR_args(parser, bed_file=True)
    args = parser.parse_args()
    return args


def main():
    """
    Convert beta file[s] to bigwig file[s].
    Assuming bedGraphToBigWig is installed and in PATH
    """
    args = parse_args()
    validate_file_list(args.beta_paths)
    if not check_executable('bedGraphToBigWig', verbose=True):
        return

    b = BetaToBigWig(args)
    for beta in args.beta_paths:
        b.run_beta_to_bw(beta)


if __name__ == '__main__':
    main()

