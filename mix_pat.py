#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_files_list, IllegalArgumentError, splitextgz, add_GR_args, delete_or_skip, BedFileWrap
from genomic_region import GenomicRegion
from merge import merge_sort_open_pats, MergePats
from pat2beta import pat2beta
from beta_cov import beta_cov, beta_cov_by_bed
import numpy as np
import sys
import pandas as pd
from view import ViewPat
import os.path as op
import subprocess


class Mixer:

    def __init__(self, pats, rates, dest_cov, gr, outdir, bed, force, strict, fast, reps, prefix):
        self.gr = gr
        self.pats = pats
        self.fast = fast
        self.reps = reps    # todo: implement reps
        self.strict = strict
        self.dest_cov = dest_cov
        self.bed_path = bed
        self.bed = None if not bed else BedFileWrap(bed)
        self.outdir = outdir
        self.force = force
        self.stats = pd.DataFrame(index=[splitextgz(f)[0] for f in self.pats])
        self.nr_pats = len(pats)
        self.dest_rates = self.validate_rates(rates)

        self.covs = self.read_covs()

        self.adj_rates = self.adjust_rates()

    def generate_output_path(self, outdir, prefix):
        if prefix:
            res = prefix
        else:
            # compose output path:
            pats_bnames = [splitextgz(op.basename(f))[0] for f in self.pats]
            res = '_'.join([str(x) for t in zip(pats_bnames, self.dest_rates) for x in t])
            region = '' if self.gr.sites is None else '_{}'.format(self.gr.region_str)
            res += '_cov_{:.2f}{}.pat'.format(self.dest_cov, region)
            res = op.join(outdir, res)
        return res

    def sample_file(self, i):
        file = self.pats[i]
        rate = self.adj_rates[i]
        out_name = '{}_rate_{:.4}_temp.pat'.format(splitextgz(op.basename(file))[0], rate)
        out_name = op.join(self.outdir, out_name)
        print('sampling {} from {}'.format(rate, file), file=sys.stderr)
        with open(out_name, 'w') as f:
            vp = ViewPat(file, f, self.gr, sub_sample=rate, bed_wrapper=self.bed, strict=self.strict)
            vp.view_pat(awk=(self.gr.sites is None))

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


        merged_path_nogz = op.join(self.outdir, name)
        if not delete_or_skip(merged_path_nogz + '.gz', self.force):
            return

        if self.fast:
            return self.fast_mix(merged_path_nogz)
        else:
            return self.slow_mix(merged_path_nogz)

    def slow_mix(self, out_nogz):
        # sample from pat files to temp files:  # todo: multiprocess it?
        tmp_files = [self.sample_file(i) for i in range(self.nr_pats)]

        self.reads_stats(tmp_files)

        # merge-sort the sampled files:
        merge_sort_open_pats(tmp_files, out_nogz, remove_open_pats=True)

    def fast_mix(self, out_nogz):
        view_flags = []
        for i in range(self.nr_pats):
            v = ' --awk '
            if self.strict:
                v += ' --strict'
            if self.bed_path is not None:
                v += ' -L {}'.format(self.bed_path)
            elif self.gr.sites is not None:
                v += ' -s {}-{}'.format(*self.gr.sites)
            v += ' --sub_sample {}'.format(self.adj_rates[i])
            view_flags.append(v)
        m = MergePats(self.pats, out_nogz)
        m.fast_merge_pats(view_flags=view_flags, label=True)

    def validate_rates(self, rates):
        if len(rates) == self.nr_pats - 1:
            rates.append(1.0 - np.sum(rates))

        if len(rates) != self.nr_pats:
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
        for i in range(self.nr_pats):
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

    parser.add_argument('--reps', type=int, default=1)

    parser.add_argument('--rates', type=float, metavar='(0.0, 1.0)', nargs='+', required=True,  #todo allow 0 and/or 1
                        help='Rates for each of the pat files. Note: the order matters!'
                             'Rate of for the last file may be omitted. '
                             'The rates will be adjusted s.t the output will be of the requested coverage.')

    parser.add_argument('--labels', nargs='+',
                        help='labels for')

    parser.add_argument('-L', '--bed_file',
                        help='Only output reads overlapping the input BED FILE. ')
    parser.add_argument('--fast', action='store_true', help='Use bash for better performance')
    parser.add_argument('--strict', action='store_true', help='Truncate reads that start/end outside the given region. '
                                                              'Only relevant if "region", "sites" '
                                                              'or "bed_file" flags are given.')

    out_or_pref = parser.add_mutually_exclusive_group(help='Provide prefix or output directory, but not both.'
                                                           'If none provided, default is current directory')
    out_or_pref.add_argument('-p', '--prefix', help='Prefix of output file')
    out_or_pref.add_argument('-o', '--out_dir', help='Output directory [.]', default='.')

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
              force=args.force,
              strict=args.strict,
              fast = args.fast,
              reps=args.reps,
              prefix=args.prefix)
    m.mix()
    m.print_rates()
    return


if __name__ == '__main__':
    main()
