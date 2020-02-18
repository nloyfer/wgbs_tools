#!/usr/bin/python3 -u

import argparse
import os.path as op
import subprocess
from utils_wgbs import delete_or_skip, load_beta_data, validate_files_list, load_dict, GenomeRefPaths, beta2vec, eprint
import os

BG_EXT = '.bedGraph'
BW_EXT = '.bigwig'
COV_BG_EXT = '_cov' + BG_EXT
COV_BW_EXT = '_cov' + BW_EXT
DEBUG_NR = 1000


class BetaToBigWig:
    def __init__(self, args):
        self.args = args
        self.debug = args.debug
        self.outdir = args.outdir
        self.chrom_sizes = GenomeRefPaths(args.genome).chrom_sizes
        self.ref_dict = self.load_dict()

    def load_dict(self):
        eprint('loading dict...')
        nrows = DEBUG_NR if self.debug else None
        rf = load_dict(nrows=nrows, genome_name=self.args.genome)
        rf['end'] = rf['start'] + 1
        rf['start'] = rf['start'] - 1
        return rf

    def bed_graph_to_bigwig(self, bed_graph, bigwig):

        # Convert bedGraph to bigWig:
        subprocess.check_call(['bedGraphToBigWig', bed_graph, self.chrom_sizes, bigwig])

        # compress or delete the bedGraph:
        if self.args.bedGraph:
            subprocess.check_call(['gzip', bed_graph])
        else:
            os.remove(bed_graph)

    def run_beta_to_bed(self, beta_path):
        eprint('{}'.format(op.basename(beta_path)))
        prefix = op.join(self.outdir, op.splitext(op.basename(beta_path))[0])
        out_bed = prefix + '.bed'
        if delete_or_skip(out_bed, self.args.force):
            return

        barr = load_beta_data(beta_path, sites=(1, DEBUG_NR + 1) if self.debug else None)
        assert (barr.shape[0] == self.ref_dict.shape[0])

        # paste dict with beta, then dump
        self.ref_dict['meth'] = barr[:, 0]
        self.ref_dict['total'] = barr[:, 1]
        self.ref_dict[self.ref_dict['total'] > 0].to_csv(out_bed, sep='\t', header=None, index=None)
        del self.ref_dict['meth'], self.ref_dict['total']

    def run_beta_to_bw(self, beta_path):
        eprint('{}'.format(op.basename(beta_path)))

        # Check if the current file should be skipped:
        prefix = op.join(self.outdir, op.splitext(op.basename(beta_path))[0])
        out_bigwig = prefix + BW_EXT
        out_bed_graph = prefix + BG_EXT
        cov_bigwig = prefix + COV_BW_EXT
        cov_bed_graph = prefix + COV_BG_EXT
        if not delete_or_skip(out_bigwig, self.args.force):
            return

        if not op.isdir(self.outdir):
            eprint('Invalid output directory:', self.outdir)
            return

        # load beta file
        barr = load_beta_data(beta_path, sites=(1, DEBUG_NR + 1) if self.debug else None)
        assert (barr.shape[0] == self.ref_dict.shape[0])

        # dump coverage:
        if self.args.dump_cov:
            eprint('Dumping cov...')
            self.ref_dict['cov'] = barr[:, 1]
            sort_and_dump_df(self.ref_dict[self.ref_dict['cov'] >= self.args.min_cov], cov_bed_graph)
            del self.ref_dict['cov']
            # convert bedGraph to bigWig:
            self.bed_graph_to_bigwig(cov_bed_graph, cov_bigwig)

        # dump beta values
        eprint('Dumping beta vals...')
        self.ref_dict['beta'] = beta2vec(barr, na=-1)
        sort_and_dump_df(self.ref_dict, out_bed_graph)
        del self.ref_dict['beta']

        # convert bedGraphs to bigWigs:
        self.bed_graph_to_bigwig(out_bed_graph, out_bigwig)


def sort_and_dump_df(df, path):
    df.sort_values(by=['chr', 'start']).to_csv(path, sep='\t', header=None, index=None)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('beta_paths', nargs='+')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-b', '--bedGraph', action='store_true', help='Keep (gzipped) bedGraphs as well as bigwigs')
    parser.add_argument('--dump_cov', action='store_true',
                        help='Generate coverage bigiwig in addition to beta values bigwig')
    parser.add_argument('-c', '--min_cov', type=int, default=1,
                        help='Minimal coverage to consider when computing beta values.'
                             ' Default is 1 (include all observations)')
    parser.add_argument('--outdir', '-o', default='.')
    parser.add_argument('--genome', help='Genome reference name. Default is hg19.', default='hg19')
    args = parser.parse_args()
    return args


def main():
    """
    Convert beta file[s] to bigwig file[s].
    """
    args = parse_args()
    validate_files_list(args.beta_paths, '.beta')

    b = BetaToBigWig(args)
    for beta in args.beta_paths:
        b.run_beta_to_bw(beta)


if __name__ == '__main__':
    main()
