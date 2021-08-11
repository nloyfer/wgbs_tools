#!/usr/bin/python3 -u

import argparse
import os
import os.path as op
import subprocess
import sys
import numpy as np
import pandas as pd

from beta_to_blocks import load_blocks_file
from betas_to_table import get_table, groups_load_wrap
from dmb import load_gfile_helper, match_prefix_to_bin, set_hypo_hyper
from utils_wgbs import validate_single_file, add_multi_thread_args, \
        validate_file_list, eprint, drop_dup_keep_order


def load_group_file(groups_file, betas):
    validate_single_file(groups_file)
    validate_file_list(betas)
    gf = load_gfile_helper(groups_file)
    gf['full_path'] = match_prefix_to_bin(gf['fname'], betas, '.beta')
    return gf


def get_validate_targets(targets, groups):
    # validate target is in groups file
    if not targets:
        return groups
    for target in targets:
        if target not in groups:
            eprint(f'Invalid target: {target}')
            eprint('Possible targets:', groups)
            raise IllegalArgumentError()
    return targets


def set_bg_tg_names(gf, targets):
    r = {}
    for group in targets:
        tg_names = gf[gf['group'] == group]['fname'].values
        bg_names = drop_dup_keep_order(gf[gf['group'] != group]['fname'])

        # remove from bg_names samples shared with tg_names
        bg_names = [s for s in bg_names if s not in tg_names]
        assert(len(bg_names) + len(tg_names) == len(set(gf['fname'])))
        r[group] = tg_names, bg_names
    return r


class MarkerFinder:
    def __init__(self, args):
        self.args = args
        self.df = pd.DataFrame()
        self.verbose = args.verbose
        self.tg_names = []
        self.bg_names = []
        self.inds_dict = {}
        self.hyper, self.hypo = set_hypo_hyper(args.hyper, args.hypo)
        self.validate_args()

        # validate output dir:
        if not op.isdir(args.out_dir):
            os.mkdir(args.out_dir)

        # load groups
        self.gf = load_group_file(args.groups_file, args.betas)
        groups = sorted(self.gf['group'].unique())
        self.targets = get_validate_targets(self.args.targets, groups)
        self.res = {t: pd.DataFrame() for t in self.targets}

    def run(self):

        # load all data
        blocks = self.load_blocks()
        self.inds_dict = set_bg_tg_names(self.gf, self.targets)
        chunk_size = self.args.chunk_size
        for start in range(0, blocks.shape[0], chunk_size):
            self.proc_chunk(blocks.iloc[start:start + chunk_size])

        # dump results:
        for target in self.targets:
            self.group = target
            self.dump_results(self.res[target].reset_index(drop=True))

    def proc_chunk(self, blocks_df):
        self.df = self.load_data_chunk(blocks_df)

        for group in self.targets:
            eprint(group)
            self.group = group
            tf = self.find_markers_group()
            self.res[group] = pd.concat([self.res[group], tf])

    def validate_args(self):
        if self.args.min_cpg < 0:
            raise IllegalArgumentError('min_cpg must be non negative')
        if self.args.max_cpg < 1:
            raise IllegalArgumentError('max_cpg must larger than 0')
        if self.args.min_bp < 0:
            raise IllegalArgumentError('min_bp must be non negative')
        if self.args.max_bp < 2:
            raise IllegalArgumentError('max_bp must larger than 1')
        validate_single_file(self.args.blocks_path)
        validate_single_file(self.args.groups_file)

    #############################
    #                           #
    #       Load methods        #
    #                           #
    #############################

    def load_blocks(self):
        # load blocks file and filter it by CpG and bg length

        df = load_blocks_file(self.args.blocks_path)
        orig_nr_blocks = df.shape[0]

        # filter by lenCpG
        df['lenCpG'] = df['endCpG'] - df['startCpG']
        df = df[df['lenCpG'] >= self.args.min_cpg]
        df = df[df['lenCpG'] <= self.args.max_cpg]

        # filter by len in bp
        df['len'] = df['end'] - df['start']
        df = df[df['len'] >= self.args.min_bp]
        df = df[df['len'] <= self.args.max_bp]

        df.reset_index(drop=True, inplace=True)

        # print stats
        if self.verbose:
            eprint(f'loaded {orig_nr_blocks:,} blocks')
            if df.shape[0] != orig_nr_blocks:
                eprint(f'droppd to {df.shape[0]:,} ')

        return df

    def load_data_chunk(self, blocks_df):
        # load methylation data from beta files collapsed to the blocks in blocks_df
        if self.verbose:
            nr_samples = len(self.gf['fname'].unique())
            eprint(f'loading data for {blocks_df.shape[0]:,} blocks over' \
                   f' {nr_samples} samples...')
        return get_table(blocks_df=blocks_df.copy(),
                         gf=self.gf,
                         min_cov=self.args.min_cov,
                         threads=self.args.threads,
                         verbose=False,
                         group=False)


    #############################
    #                           #
    #       Finding markers     #
    #                           #
    #############################

    def find_markers_group(self):
        if self.df.empty:
            return self.df

        # load context
        self.tg_names, self.bg_names = self.inds_dict[self.group]
        tf = self.df.copy()

        # filter blocks by coverage:
        df_tg = tf[self.tg_names]
        df_bg = tf[self.bg_names]
        keep_tg = (df_tg.notna().sum(axis=1) / (df_tg.shape[1])) >= self.args.na_thresh
        keep_bg = (df_bg.notna().sum(axis=1) / (df_bg.shape[1])) >= self.args.na_thresh
        tf = tf.loc[keep_tg & keep_bg, :].reset_index(drop=True)

        # find markers:
        tfM = self.find_M_markers(tf)
        tfU = self.find_U_markers(tf)
        tf = pd.concat([tfU, tfM]).reset_index(drop=True)
        tf['delta'] = ((tf['tg'] - tf['bg']).abs() * 100).astype(int)
        return tf

    def find_M_markers(self, tf):
        # look for 'M' markers
        if not self.hyper:
            return pd.DataFrame()

        df_bg = tf[self.bg_names].fillna(.5)
        keep = df_bg.mean(axis=1) <= self.args.bg_beta_M
        tfM = tf.loc[keep, :].copy()
        df_tg = tf.loc[keep, self.tg_names].fillna(.5)
        df_bg = tf.loc[keep, self.bg_names].fillna(.5)

        tfM['tg'] = np.quantile(df_tg,
                                self.args.tg_quant,
                                interpolation='higher',
                                axis=1)
        tfM['bg'] = np.quantile(df_bg,
                                1 - self.args.bg_quant,
                                interpolation='lower',
                                axis=1)
        tfM = tfM[(tfM['tg'] - tfM['bg'] >= self.args.delta)]
        tfM['direction'] = 'M'
        return tfM

    def find_U_markers(self, tf):
        # look for 'U' markers
        if not self.hypo:
            return pd.DataFrame()

        df_bg = tf[self.bg_names].fillna(.5)
        keep = df_bg.mean(axis=1) >= self.args.bg_beta_U
        tfU = tf.loc[keep, :].copy()
        df_tg = tf.loc[keep, self.tg_names].fillna(.5)
        df_bg = tf.loc[keep, self.bg_names].fillna(.5)

        tfU['tg'] = np.quantile(df_tg,
                                1 - self.args.tg_quant,
                                interpolation='lower',
                                axis=1)
        tfU['bg'] = np.quantile(df_bg,
                                self.args.bg_quant,
                                interpolation='higher',
                                axis=1)
        tfU = tfU[(tfU['bg'] - tfU['tg'] >= self.args.delta)]
        tfU['direction'] = 'U'
        return tfU



    #############################
    #                           #
    #       Dump methods        #
    #                           #
    #############################

    def dump_results(self, tf):
        eprint(f'Number of markers found: {tf.shape[0]:,}')
        if self.args.top:
            tf = tf.sort_values(by='delta', ascending=False).head(self.args.top)
        tf = tf.sort_values(by='startCpG')
        tf['target'] = self.group
        tf['lenCpG'] = tf['lenCpG'].astype(str) + 'CpGs'
        tf['bp'] = (tf['end'] - tf['start']).astype(str) + 'bp'
        tf['region'] = tf['chr'] + ':' + tf['start'].astype(str) + '-' + tf['end'].astype(str)
        cols_to_dump = ['chr', 'start', 'end', 'startCpG', 'endCpG',
                        'target', 'region', 'lenCpG', 'bp', 'tg',
                        'bg', 'delta', 'direction']
        self.dump_single_df(tf[cols_to_dump].round(2))

    def dump_single_df(self, df):
        """ Dump the results, with comments header """

        outpath = op.join(self.args.out_dir, f'Markers.{self.group}.bed')
        eprint(f'dumping to {outpath}')
        df_mode = 'w'

        # write header:
        if self.args.header:
            # write header (comments):
            with open(outpath, 'w') as f:
                for sample in sorted(self.tg_names):
                    f.write(f'#> {sample}\n' )
                for sample in sorted(self.bg_names):
                    f.write(f'#< {sample}\n' )
            df_mode = 'a'

        # dump df:
        df.to_csv(outpath, index=None, sep='\t', mode=df_mode, header=None)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--blocks_path', '-b', help='Blocks bed path.', required=True)
    parser.add_argument('--groups_file', '-g', help='csv file of groups', required=True)
    parser.add_argument('--targets', nargs='+', help='find markers only for these groups')
    parser.add_argument('--betas', help='beta file paths. files not in the group files are ignored', nargs='+', required=True)
    parser.add_argument('-o', '--out_dir', default='.', help='Output directory [.]')
    parser.add_argument('--min_bp', type=int, default=0)
    parser.add_argument('--max_bp', type=int, default=1000000)
    parser.add_argument('--min_cpg', type=int, default=0)
    parser.add_argument('--max_cpg', type=int, default=1000000)
    parser.add_argument('-m', '--delta', type=float, default=0.5,
                       help='Filter markers by beta values delta. range: [0.0, 1.0]. Default 0.5')
    parser.add_argument('-c', '--min_cov', type=int, default=5,
            help='Minimal number of binary observations in block coverage to be considered [5]')
    parser.add_argument('--hyper', action='store_true', help='Only consider hyper-methylated markers')
    parser.add_argument('--hypo', action='store_true', help='Only consider hypo-methylated markers')
    parser.add_argument('--top', type=int,
                        help='Output only the top TOP markers, under the constraints. [All]')
    parser.add_argument('--header', action='store_true', help='add header to output files')
    parser.add_argument('--tg_quant', type=float, default=.25, help='quantile of target samples to ignore')
    parser.add_argument('--bg_quant', type=float, default=.025, help='quantile of background samples to ignore')

    parser.add_argument('--bg_beta_U', type=float, default=.6, help='background average beta value for "U" markers [.6]')
    parser.add_argument('--bg_beta_M', type=float, default=.5, help='background average beta value for "M" markers [.5]')

    parser.add_argument('--na_thresh', type=float, default=.66, help='rate of samples with sufficient coverage required in both target and background [.66]')

    parser.add_argument('--chunk_size', type=int, default=150000, help='Number of blocks to load on each step [150000]')
    parser.add_argument('--verbose', '-v', action='store_true')
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args



def main():
    """
    Find differentially methylated blocks
    """
    args = parse_args()
    MarkerFinder(args).run()


if __name__ == '__main__':
    main()
