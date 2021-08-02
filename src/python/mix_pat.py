#!/usr/bin/python3 -u

import argparse
import numpy as np
import pandas as pd
import os.path as op
from multiprocessing import Pool
import multiprocessing
from src.python.utils_wgbs import validate_file_list, IllegalArgumentError, splitextgz, add_GR_args, delete_or_skip, \
        BedFileWrap, eprint, validate_dir, add_multi_thread_args
from src.python.genomic_region import GenomicRegion
from src.python.merge import MergePats
from src.python.pat2beta import pat2beta
from src.python.beta_cov import beta_cov, beta_cov_by_bed


class Mixer:

    def __init__(self, args):
        eprint('mixing...')
        self.args = args
        self.gr = GenomicRegion(args)
        self.pats = args.pat_files
        self.dest_cov = args.cov
        self.bed = None if not args.bed_file else BedFileWrap(args.bed_file)
        self.stats = pd.DataFrame(index=[splitextgz(op.basename(f))[0] for f in self.pats])
        self.nr_pats = len(self.pats)
        self.labels = self.validate_labels(args.labels)

        self.dest_rates = self.validate_rates(args.rates)
        self.covs = self.read_covs()
        self.adj_rates = self.adjust_rates()

        self.prefix = self.generate_prefix(args.out_dir, args.prefix)

    def generate_prefix(self, outdir, prefix):
        if prefix:
            if op.dirname(prefix):
                validate_dir(op.dirname(prefix))
            return prefix
        else:
            validate_dir(outdir)
            # compose output path:
            pats_bnames = [splitextgz(op.basename(f))[0] for f in self.pats]
            res = '_'.join([str(x) for t in zip(pats_bnames, self.dest_rates) for x in t])
            region = '' if self.gr.sites is None else '_{}'.format(self.gr.region_str)
            res += '_cov_{:.2f}{}'.format(self.dest_cov, region)
            res = op.join(outdir, res)
        return res

    def print_rates(self):
        eprint('Requested Coverage: {:.2f}'.format(self.dest_cov))
        eprint(self.stats)

    def add_stats_col(self, title, data):
        self.stats[title] = data

    def single_mix(self, rep):
        prefix_i = self.prefix + '_{}.pat'.format(rep + 1)
        if not delete_or_skip(prefix_i + '.gz', self.args.force):
            return

        view_flags = []
        for i in range(self.nr_pats):
            # v = ' --awk '
            v = ''
            if self.args.strict:
                v += ' --strict'
            if self.args.bed_file is not None:
                v += ' -L {}'.format(self.args.bed_file)
            elif self.gr.sites is not None:
                v += ' -s {}-{}'.format(*self.gr.sites)
            v += ' --sub_sample {}'.format(self.adj_rates[i])
            view_flags.append(v)
        eprint('prefix:', prefix_i)
        m = MergePats(self.pats, prefix_i, self.labels, args=self.args)
        m.fast_merge_pats(view_flags=view_flags)

    def validate_labels(self, labels):
        if labels is None:
            labels = [splitextgz(op.basename(p))[0].split('-')[0].lower() for p in self.pats]

        if len(labels) != self.nr_pats:
            raise IllegalArgumentError('len(labels) != len(files)')
        return labels

    def validate_rates(self, rates):
        if len(rates) == self.nr_pats - 1:
            rates.append(1.0 - np.sum(rates))

        if len(rates) != self.nr_pats:
            raise IllegalArgumentError('len(rates) must be in {len(files), len(files) - 1}')

        if np.abs(np.sum(rates) - 1) > 1e-8:
            raise IllegalArgumentError('Sum(rates) == {} != 1'.format(np.sum(rates)))

        if np.min(rates) < 0 or np.max(rates) > 1:
            raise IllegalArgumentError('rates must be in range [0, 1)')

        self.add_stats_col('ReqstRates', rates)
        return rates

    def read_covs(self):
        covs = []
        for pat in self.pats:
            suff = '.lbeta' if self.args.lbeta else '.beta'
            beta = pat.replace('.pat.gz', suff)
            if not op.isfile(beta):
                eprint('No {} file compatible to {} was found. Generate it...'.format(suff, pat))
                pat2beta(pat, op.dirname(pat), args=self.args, force=True)
            if self.bed:
                cov = beta_cov_by_bed(beta, self.bed)
            elif self.args.bed_cov:     # todo: this is messy. fix it
                cov = beta_cov_by_bed(beta, BedFileWrap(self.args.bed_cov))
            else:
                cov = beta_cov(beta, self.gr.sites, print_res=True)
            covs.append(cov)
        self.add_stats_col('OrigCov', covs)
        return covs

    def adjust_rates(self):

        if not self.dest_cov:
            self.dest_cov = self.covs[int(np.argmax(self.dest_rates))]

        adj_rates = []
        for i in range(self.nr_pats):
            adjr = self.dest_rates[i] * self.dest_cov / self.covs[i]
            if adjr > 1:
                self.print_rates()
                raise IllegalArgumentError('File {} has too low coverage'.format(self.pats[i]))
            adj_rates.append(adjr)

        self.add_stats_col('AdjRates', adj_rates)
        return adj_rates


def single_mix(i, m):
    m.single_mix(i)


def mult_mix(args):
    m = Mixer(args)
    m.print_rates()
    with Pool(args.threads) as p:
        for i in range(args.reps):
            p.apply_async(single_mix, (i, m))
        p.close()
        p.join()
    # m.print_rates()


##########################
#                        #
#         Main           #
#                        #
##########################

# todo: support unq
def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat_files', nargs='+', help='Two or more pat files')
    parser.add_argument('--bed_cov', help='calculate coverage on this bed file regions only') # todo: remove or validate file exists etc.
    parser.add_argument('-c', '--cov', type=float,
                        help='Coverage of the output pat. '
                             'Default the coverage of the file with the highest rate. '
                             'Only supported if corresponding beta files are in the same '
                             'directory with the pat files. '
                             'Otherwise, they will be created.')
    add_GR_args(parser, bed_file=True)
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')

    parser.add_argument('--reps', type=int, default=1, help='nr or repetitions [1]')

    parser.add_argument('--rates', type=float, metavar='[0.0, 1.0]', nargs='+', required=True,
                        help='Rates for each of the pat files. Note: the order matters!'
                             'Rate of for the last file may be omitted. '
                             'The rates will be adjusted s.t the output will be of the requested coverage.')

    parser.add_argument('--labels', nargs='+', help='labels for the mixed reads. '
                                                    'Default is the basenames of the pat files,'
                                                    'lowercased and trimmed by the first "-"')

    parser.add_argument('--strict', action='store_true', help='Truncate reads that start/end outside the given region. '
                                                              'Only relevant if "region", "sites" '
                                                              'or "bed_file" flags are given.')

    out_or_pref = parser.add_mutually_exclusive_group()
    out_or_pref.add_argument('-p', '--prefix', help='Prefix of output file.')
    out_or_pref.add_argument('-o', '--out_dir', help='Output directory [.]', default='.')
    parser.add_argument('-T', '--temp_dir', help='passed to "sort -m". Useful for merging very large pat files')
    parser.add_argument('-l', '--lbeta', action='store_true', help='Use lbeta file (uint16) instead of beta (uint8)')
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args


def main():
    """
    Mix samples from K different pat files.
    Output a single mixed pat.gz[.csi] file - sorted, bgzipped and indexed -
    with an informative name.
    """
    args = parse_args()
    validate_file_list(args.pat_files, 'pat.gz', 2)
    mult_mix(args)
    return


if __name__ == '__main__':
    main()
