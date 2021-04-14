#!/usr/bin/python3 -u

import argparse
import subprocess
import numpy as np
import sys
import os.path as op
import pandas as pd
import math
import os
from tqdm import tqdm
from utils_wgbs import load_beta_data2, validate_single_file, eprint, \
                       IllegalArgumentError, GenomeRefPaths

DEFAULT_BLOCKS = GenomeRefPaths().blocks
DEBUG_NR = 100000
um_ind_dict = {'U': 0, 'M': 2}


class MarkersFinder:
    def __init__(self, args):
        self.args = args
        self.dfU = pd.DataFrame()
        self.dfM = pd.DataFrame()
        self.blocks = pd.DataFrame()
        self.nr_blocks = 0
        self.orig_nr_blocks = 0
        self.keepinds = None
        self.groups = None
        self.verbose = args.verbose
        self.hyper, self.hypo = self.set_hypo_hyper()
        self.validate_args()
        # validate output dir:
        if not op.isdir(args.out_dir):
            os.mkdir(args.out_dir)
        # load groups
        self.gf = load_groups_file(args.groups_file, args.input_dir, args.verbose)
        self.gf_nodup = self.gf.drop_duplicates(subset='fname').reset_index(drop=True)
        # validate target is in groups file
        target = self.args.target
        if target and target not in self.gf['group'].values:
            eprint(f'target {target} not in groups file {self.args.groups_file}')
            eprint('Possible targets:', sorted(self.gf['group'].unique()))
            raise IllegalArgumentError()

    def set_hypo_hyper(self):
        hyper, hypo = self.args.hyper, self.args.hypo
        if not hyper and not hypo:
            return True, True
        return hyper, hypo

    def run(self):

        # load all data
        self.blocks = self.load_blocks_file()
        self.dfU, self.dfM = self.load_bins()
        for group in sorted(self.gf['group'].unique()):
            if self.args.target and group != self.args.target:
                continue
            eprint(group)
            self.group = group
            tfU = self.find_markers_group('U')
            tfM = self.find_markers_group('M')
            tf = pd.concat([tfU, tfM])
            self.dump_results(tf)

    def validate_args(self):
        if self.args.min_cpg < 1:
            raise IllegalArgumentError('min_cpg must be a positive integer')
        validate_single_file(self.args.blocks_path)
        validate_single_file(self.args.groups_file)

    #############################
    #                           #
    #       Load methods        #
    #                           #
    #############################

    def load_blocks_file(self):
        if self.verbose:
            eprint('loading blocks...')
        names = ['chr', 'start', 'end', 'startCpG', 'endCpG']
        cols = range(len(names))
        nrows = DEBUG_NR if self.args.debug else None
        df = pd.read_csv(self.args.blocks_path, sep='\t',
                header=None, names=names, nrows=nrows, usecols=cols)
        df['lenCpG'] = df['endCpG'] - df['startCpG']
        self.keepinds = df['lenCpG'] >= self.args.min_cpg
        self.orig_nr_blocks = df.shape[0]
        df = df[self.keepinds].reset_index(drop=True)
        self.nr_blocks = df.shape[0]
        if self.verbose:
            eprint(f'loaded {self.orig_nr_blocks:,}')
            if self.nr_blocks != self.orig_nr_blocks:
                eprint(f'droppd to {self.nr_blocks:,} with >={self.args.min_cpg} CpGs')
        return df

    def load_bins(self):
        if self.verbose:
            eprint('loading bins...')
        # breakpoint()
        nr_cols = (3 if self.args.uxm else 2)
        binsize = self.gf['binsize'][0] / self.orig_nr_blocks
        binsize /= nr_cols
        if binsize != int(binsize):
            raise IllegalArgumentError('Error: bin file size does not match blocks number')

        dtype = np.uint8 if binsize == 1 else np.uint16

        dfU = pd.DataFrame()
        dfM = pd.DataFrame()
        if self.hypo:
            dfU = np.zeros((self.nr_blocks, self.gf_nodup.shape[0]), dtype=np.float16)
        if self.hyper:
            dfM = np.zeros((self.nr_blocks, self.gf_nodup.shape[0]), dtype=np.float16)

        for ind, row in tqdm(self.gf_nodup.iterrows(), total=self.gf_nodup.shape[0]):
            data = np.fromfile(row['full_path'], dtype).reshape((-1, nr_cols))[self.keepinds, :]
            if self.hypo:
                dfU[:, ind] = self.table2vec(data, 'U')
            if self.hyper:
                dfM[:, ind] = self.table2vec(data, 'M')

        return self.array2df(dfU), self.array2df(dfM)

    def array2df(self, df):
        if not df.shape[0]:
            return pd.DataFrame()
        return pd.concat([self.blocks,
                          pd.DataFrame(data=df, columns=self.gf_nodup['fname'].values)],
                          axis=1)

    def table2vec(self, data, um):
        covs = data.sum(axis=1)
        cond = covs > self.args.min_cov
        r = np.divide(data[:, um_ind_dict[um]], covs, where=cond)
        r[~cond] = np.nan
        return r.astype(np.float16)

    def find_markers_group(self, um):
        df = self.dfU if um == 'U' else self.dfM
        if df.empty:
            return pd.DataFrame(columns=self.blocks.columns)
        self.tg_names = self.gf[self.gf['group'] == self.group]['fname'].values
        bg_names = self.gf[self.gf['group'] != self.group]['fname'].unique()

        # remove from bg_names samples shared with tg_names
        self.bg_names = [s for s in bg_names if s not in self.tg_names]
        assert(len(self.bg_names) + len(self.tg_names) == self.gf_nodup.shape[0])

        tf = self.blocks.copy()

        tf[f'tg_low'] = np.quantile(df[self.tg_names].fillna(.5),
                self.args.tg_quant, interpolation='higher', axis=1)
        tf[f'bg_high'] = np.quantile(df[self.bg_names].fillna(.5),
                1 - self.args.bg_quant, interpolation='lower', axis=1)
        tf = tf[(tf['tg_low'] - tf['bg_high'] >= self.args.margin)]
        tf['margin'] = ((tf['tg_low'] - tf['bg_high']) * 100).astype(int)
        tf['direction'] = um
        return tf


    #############################
    #                           #
    #       Dump methods        #
    #                           #
    #############################

    def dump_results(self, tf):
        eprint(f'Number of markers found: {tf.shape[0]:,}')
        if self.args.top:
            tf = tf.head(self.args.top)
        tf = tf.sort_values(by='startCpG', ascending=True)
        tf['target'] = self.group
        tf['lenCpG'] = tf['lenCpG'].astype(str) + 'CpGs'
        tf['bp'] = (tf['end'] - tf['start']).astype(str) + 'bp'
        tf['region'] = tf['chr'] + ':' + tf['start'].astype(str) + '-' + tf['end'].astype(str)
        tf = tf.round(2)
        outpath = op.join(self.args.out_dir, f'Markers.{self.group}.bed')
        cols_to_dump = ['chr', 'start', 'end', 'startCpG', 'endCpG',
                        'target', 'region', 'lenCpG', 'bp', 'tg_low',
                        'bg_high', 'direction', 'margin']
        self.dump_single_df(tf[cols_to_dump], outpath)

    def dump_single_df(self, df, outpath):
        """ Dump the results, with comments header """

        eprint(f'dumping to {outpath}')
        df_mode = 'w'
        # write header:
        if not self.args.no_header:
            # write header (comments):
            with open(outpath, 'w') as f:
                for sample in sorted(self.tg_names):
                    f.write(f'#> {sample}\n' )
                for sample in sorted(self.bg_names):
                    f.write(f'#< {sample}\n' )
            df_mode = 'a'

        # dump df:
        df.to_csv(outpath, index=None, sep='\t', mode=df_mode, header=None)

def load_gfile_helper(groups_file):
    # load and validate csv
    gf = pd.read_csv(groups_file, index_col=False, comment='#')
    if 'group' not in gf.columns:
        raise IllegalArgumentError('gropus file must have a column named "group"')
    # drop samples where include==False
    if 'include' in gf.columns:
        gf = gf[gf['include']]
    # drop unnecessary columns
    gf = gf.rename(columns={gf.columns[0]: 'fname'})
    gf = gf[['fname', 'group']].dropna().reset_index(drop=True)
    return gf

def load_groups_file(groups_file, idir, verbose=False):
    if verbose:
        eprint('loading groups_file...')
    gf = load_gfile_helper(groups_file)
    # find binary file paths
    return find_bin_paths(gf, idir)

def match_prefix_to_bin(prefixes, bins):
    full_paths = []
    for prefix in prefixes:
        # look for bin files starting with current prefix:
        results = [f for f in bins if op.basename(f).startswith(prefix)]
        # make sure there is exactly one such bin file
        if not results:
            raise IllegalArgumentError(f'Invalid prefix: {prefix}. Not in input bins')
        if len(results) > 1:
            raise IllegalArgumentError(f'Found multiple matches for prefix {prefix}')
        full_paths.append(results[0])
    return full_paths

def find_bin_paths(gf, idir):
    # find full path for all samples
    if not op.isdir(idir):
        raise IllegalArgumentError(f'Invalid input directory: {idir}')
    all_bins = [f for f in os.listdir(idir) if f.endswith(('.bin', '.lbeta', '.uxm', '.beta'))]

    gf['full_path'] = [op.join(indir, f) for f in match_prefix_to_bin(gf['fname'], all_bins)]

    # validation: make sure all bin files have the same size
    gf['binsize'] = gf['full_path'].apply(op.getsize)
    if len(gf['binsize'].unique()) != 1:
        eprint(f'Error. All binary files must have similar size')
        eprint(gf.sort_values(by='binsize', ignore_index=True).iloc[[0, -1], [0, -1]])
        raise IllegalArgumentError(f'Error. binary file sizes mismatch')

    return gf



def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--blocks_path', '-b', default=DEFAULT_BLOCKS,
            help=f'Blocks bed path. Default [{DEFAULT_BLOCKS}]')

    parser.add_argument('--groups_file', '-g', help='csv file of groups', required=True)
    parser.add_argument('--input_dir', '-i', help='directory with binary files.', required=True)

    parser.add_argument('--uxm', action='store_true')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('-o', '--out_dir', default='.', help='Output directory [.]')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-m', '--margin', type=float, default=0.5,
                       help='Filter markers by beta values margin. range: [0., 1.]')
    parser.add_argument('--tg_quant', type=float, default=.25, help='quantile of target samples to ignore')
    parser.add_argument('--bg_quant', type=float, default=.025, help='quantile of background samples to ignore')
    parser.add_argument('-c', '--min_cov', type=int, default=4, help='Minimal coverage to be considered,'
                                                                 'In both groups. [4]')
    parser.add_argument('--hyper', action='store_true', help='Only consider markers in which target group'
                                                             ' is hyper methylated')
    parser.add_argument('--hypo', action='store_true', help='Only consider markers in which target group'
                                                            ' is hypo methylated')
    parser.add_argument('--min_cpg', type=int, default=3, help='Drop blocks with less than MIN_CPG sites [3]')
    parser.add_argument('--target', help='only find markers for this group [all]')
    parser.add_argument('--top', type=int,
                        help='Output only the top TOP markers, under the constraints. [All]')
    parser.add_argument('--debug', '-d', action='store_true',
                        help=f'Debug mode. Only consider first {DEBUG_NR} blocks')
    parser.add_argument('--no_header', '-nh', action='store_true',
                        help='Don\'t output comments in the beginning of the output file ')
    args = parser.parse_args()
    return args


def main():
    """
    Find markers (blocks) to differentiate between two or more groups of samples
    (collapsed beta files or homog binary files).
    """
    args = parse_args()
    if not args.uxm:
        eprint('Only --uxm mode is currently supported')
        raise NotImplementedError    #todo: implement
    MarkersFinder(args).run()


if __name__ == '__main__':
    main()

