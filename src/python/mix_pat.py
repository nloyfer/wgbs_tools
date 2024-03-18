#!/usr/bin/python3 -u

import argparse
import os.path as op
from multiprocessing import Pool
import numpy as np
import pandas as pd
from utils_wgbs import validate_file_list, IllegalArgumentError, delete_or_skip, \
        eprint, validate_dir, add_multi_thread_args, pretty_name
from genomic_region import GenomicRegion
from merge import MergePats, validate_labels, extract_view_flags
from cview import add_view_flags
from pat2beta import pat2beta
from beta_cov import beta_cov
from beta_to_blocks import load_blocks_file


class Mixer:

    def __init__(self, args):
        self.args = args
        self.gr = GenomicRegion(args)
        self.pats = args.pat_files
        self.dest_cov = args.cov
        self.bed = load_blocks_file(args.bed_file) if args.bed_file else None
        self.stats = pd.DataFrame(index=[pretty_name(f) for f in self.pats])
        self.nr_pats = len(self.pats)

        self.dest_rates = self.validate_rates(args.rates)
        self.covs = self.read_covs()
        self.adj_rates = self.adjust_rates()

        self.prefix = self.generate_prefix(args.out_dir, args.prefix)

    def generate_prefix(self, outdir, prefix):
        if prefix:
            if op.dirname(prefix):
                validate_dir(op.dirname(prefix))
            return prefix

        validate_dir(outdir)
        # compose output path:

        pats_bnames = [pretty_name(f) for f in self.pats]
        res = '_'.join([str(x) for t in zip(pats_bnames, self.dest_rates) for x in t])
        region = '' if self.gr.sites is None else f'_{self.gr.region_str}'
        res += '_cov_{:.2f}{}'.format(self.dest_cov, region)
        res = op.join(outdir, res)
        return res

    def print_rates(self):
        eprint('Requested Coverage: {:.2f}'.format(self.dest_cov))
        eprint(self.stats)

    def add_stats_col(self, title, data):
        self.stats[title] = data

    def single_mix(self, rep):
        mix_i = self.prefix + f'_{rep + 1}.pat.gz'
        if not delete_or_skip(mix_i, self.args.force):
            return
        eprint(f'[wt mix] mix: {mix_i}')

        v = extract_view_flags(self.args)
        view_flags = [v + f' --sub_sample {s}' for s in self.adj_rates]

        m = MergePats(self.pats, mix_i,
                      validate_labels(self.args.labels, self.pats, required=True),
                      args=self.args)
        m.fast_merge_pats(view_flags=view_flags)

    def validate_rates(self, rates):
        if len(rates) == self.nr_pats - 1:
            rates.append(1.0 - np.sum(rates))

        if len(rates) != self.nr_pats:
            raise IllegalArgumentError('len(rates) must be in {len(files), len(files) - 1}')

        if np.abs(np.sum(rates) - 1) > 1e-8:
            raise IllegalArgumentError(f'Sum(rates) == {np.sum(rates)} != 1')

        if np.min(rates) < 0 or np.max(rates) > 1:
            raise IllegalArgumentError('rates must be in range [0, 1)')

        self.add_stats_col('ReqstRates', rates)
        return rates

    def read_covs(self):
        suff = '.lbeta' if self.args.lbeta else '.beta'
        covs = []
        for pat in self.pats:
            beta = pat[:-7] + suff
            if not op.isfile(beta):
                eprint(f'[wt mix] No {suff} file compatible to {pat} was found. Generating it...')
                pat2beta(pat, op.dirname(pat), args=self.args, force=True)
            cov = beta_cov(beta, blocks_df=self.bed, sites=self.gr.sites)
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
                eprint(f'[wt mix] WARNING: {self.pats[i]} has low coverage. Reads will be duplicated')
            adj_rates.append(adjr)

        self.add_stats_col('AdjRates', adj_rates)
        return adj_rates


def single_mix(i, m):
    m.single_mix(i)


def mult_mix(args):
    m = Mixer(args)
    m.print_rates()
    p = Pool(args.threads)
    params = [(i, m) for i in range(args.reps)]
    p.starmap(single_mix, params)
    p.close()
    p.join()


##########################
#                        #
#         Main           #
#                        #
##########################

def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat_files', nargs='+', help='Two or more pat files')
    # parser.add_argument('--bed_cov', help='calculate coverage on this bed ' \
                        # 'file regions only') # todo: remove or validate file exists etc.
    parser.add_argument('-c', '--cov', type=float,
                        help='Coverage of the output pat. '
                             'Default the coverage of the file with the highest rate. '
                             'Only supported if corresponding beta files are in the same '
                             'directory with the pat files. '
                             'Otherwise, they will be created.')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Overwrite existing files if existed')

    parser.add_argument('--reps', type=int, default=1, help='nr or repetitions [1]')

    parser.add_argument('--rates', type=float, metavar='[0.0, 1.0]', nargs='+', required=True,
                        help='Rates for each of the pat files. Note: the order matters!'
                             'Rate of for the last file may be omitted. '
                             'The rates will be adjusted s.t the output will be of the requested coverage.')

    parser.add_argument('--labels', nargs='+', help='labels for the mixed reads. '
                                                    'Default is the basenames of the pat files')

    out_or_pref = parser.add_mutually_exclusive_group()
    out_or_pref.add_argument('-p', '--prefix', help='Prefix of output file.')
    out_or_pref.add_argument('-o', '--out_dir', help='Output directory [.]', default='.')
    parser.add_argument('-T', '--temp_dir', help='passed to "sort -m". Useful for merging very large pat files')
    parser.add_argument('-l', '--lbeta', action='store_true', help='Use lbeta file (uint16) instead of beta (uint8)')
    parser.add_argument('-v', '--verbose', action='store_true')
    add_view_flags(parser, sub_sample=False, out_path=False)
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


if __name__ == '__main__':
    main()
