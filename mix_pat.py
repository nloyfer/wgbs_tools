#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_files_list, IllegalArgumentError, splitextgz, add_GR_args, delete_or_skip, BedFileWrap
from genomic_region import GenomicRegion
from merge import merge_sort_open_pats
from pat2beta import pat2beta
from beta_cov import beta_cov, beta_cov_by_bed
import numpy as np
import sys
import pandas as pd
from view import ViewPat
import os.path as op
import subprocess



class Mixer:

    def __init__(self, pats, rates, dest_cov, gr, outdir, bed, force):
        self.gr = gr
        self.pats = pats
        self.dest_cov = dest_cov
        self.bed = None if not bed else BedFileWrap(bed)
        self.outdir = outdir
        self.force = force
        self.stats = pd.DataFrame(index=[splitextgz(f)[0] for f in self.pats])
        self.k = len(pats)
        self.dest_rates = self.validate_rates(rates)

        self.covs = self.read_covs()

        self.adj_rates = self.adjust_rates()

    def sample_file(self, i):
        file = self.pats[i]
        rate = self.adj_rates[i]
        out_name = '{}_rate_{:.4}_temp.pat'.format(splitextgz(op.basename(file))[0], rate)
        out_name = op.join(self.outdir, out_name)
        with open(out_name, 'w') as f:
            ViewPat(file, f, self.gr, sub_sample=rate, bed_wrapper=self.bed).view_pat()
        return out_name

    def reads_stats(self, tmp_files):
        nr_lines = []
        for pat in tmp_files:
            if not op.getsize(pat):
                n = 0
            else:
                cmd = "awk '{s+=$4}END{print s}' " + pat
                n = int(subprocess.check_output(cmd, shell=True))
            nr_lines.append(n)
        nr_lines = np.array(nr_lines)
        actual_rates = nr_lines / np.sum(nr_lines)
        self.add_stats_col('ResRate', actual_rates)
        return actual_rates

    def print_rates(self):
        print('Requested Coverage: {:.2f}'.format(self.dest_cov), file=sys.stderr)
        print(self.stats, file=sys.stderr)

    def add_stats_col(self, title, data):
        self.stats[title] = data

    def mix(self):

        # compose output path:
        name = '_'.join([str(x) for t in zip([splitextgz(op.basename(f))[0] for f in self.pats], self.dest_rates) for x in t])
        r = '' if self.gr.sites is None else '_{}'.format(self.gr.region_str)
        name += '_cov_{:.2f}{}.pat'.format(self.dest_cov, r)
        merged_path_nogz = op.join(self.outdir, name)
        if not delete_or_skip(merged_path_nogz + '.gz', self.force):
            return

        # sample from pat files to temp files:  # todo: multiprocess it
        tmp_files = [self.sample_file(i) for i in range(self.k)]

        self.reads_stats(tmp_files)

        # merge-sort the sampled files:
        merge_sort_open_pats(tmp_files, merged_path_nogz, remove_open_pats=False)

    def validate_rates(self, rates):
        if len(rates) == self.k - 1:
            rates.append(1.0 - np.sum(rates))
        if len(rates) != self.k:
            raise IllegalArgumentError('len(rates) must be in {len(files), len(files) - 1}')
        if np.sum(rates) != 1:
            raise IllegalArgumentError('Sum(rates) == {} != 1'.format(np.sum(rates)))
        if np.min(rates) <= 0 or np.max(rates) >= 1:
            raise IllegalArgumentError('rates must be in range (0, 1)')
        self.add_stats_col('ReqstRates', rates)
        return rates

    def read_covs(self):
        covs = []
        for pat in self.pats:
            beta = pat.replace('.pat.gz', '.beta')
            if not op.isfile(beta):
                print('No beta file compatible to {} was found. '
                      'Generate it...'.format(pat), file=sys.stderr)
                pat2beta(pat, op.dirname(pat))
            if self.bed:
                cov = beta_cov_by_bed(beta, self.bed)
            else:
                cov = beta_cov(beta, self.gr.sites)
            covs.append(cov)
        self.add_stats_col('OrigCov', covs)
        return covs

    def adjust_rates(self):

        if not self.dest_cov:
            self.dest_cov = self.covs[int(np.argmax(self.dest_rates))]

        adj_rates = []
        for i in range(self.k):
            adjr = self.dest_rates[i] * self.dest_cov / self.covs[i]
            if adjr > 1:  # todo: allow adj_rate > 1?
                self.print_rates()
                raise IllegalArgumentError('File {} has too low coverage'.format(self.pats[i]))
            adj_rates.append(adjr)

        self.add_stats_col('AdjRates', adj_rates)
        return adj_rates


##########################
#                        #
#         Main           #
#                        #
##########################

# todo: support unq
def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat_files', nargs='+', help='Two or more pat files')
    parser.add_argument('-c', '--cov', type=float,
                        help='Coverage of the output pat. '
                             'Default the coverage of the file with the highest rate. '
                             'Only supported if corresponding beta files are in the same '
                             'directory with the pat files. '
                             'Otherwise, they will be created.')
    add_GR_args(parser)
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')
    parser.add_argument('-o', '--out_dir', help='Output directory [.]', default='.')
    parser.add_argument('--rates', type=float, metavar='(0.0, 1.0)', nargs='+', required=True,
                        help='Rates for each of the pat files. Note: the order matters!'
                             'Rate of for the last file may be omitted. '
                             'The rates will be adjusted s.t the output will be of the requested coverage.')
    parser.add_argument('-L', '--bed_file',
                        help='Only output reads overlapping the input BED FILE. ')
    args = parser.parse_args()
    return args


def main():
    """
    Mix samples from K different pat files.
    Output a single mixed pat.gz[.csi] file - sorted, bgzipped and indexed -
    with an informative name.
    """
    args = parse_args()
    validate_files_list(args.pat_files, 'pat.gz', 2)

    if args.bed_file and (args.region or args.sites):
        print('-L, -s and -r are mutually exclusive', file=sys.stderr)
        return

    m = Mixer(pats=args.pat_files,
              rates=args.rates,
              gr=GenomicRegion(args),
              dest_cov=args.cov,
              outdir=args.out_dir,
              bed=args.bed_file,
              force=args.force)
    m.mix()
    m.print_rates()
    return


if __name__ == '__main__':
    main()
