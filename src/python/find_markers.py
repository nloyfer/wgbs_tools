#!/usr/bin/python3 -u

import argparse
import subprocess
import numpy as np
import sys
import os.path as op
import pandas as pd
import math
import os
from utils_wgbs import validate_single_file, add_multi_thread_args, \
        validate_file_list, eprint, drop_dup_keep_order
from betas_to_table import get_table, groups_load_wrap
from dmb import load_gfile_helper, match_prefix_to_bin, set_hypo_hyper


DEBUG_NR = 100000


def load_group_file(groups_file, betas):
    validate_single_file(groups_file)
    validate_file_list(betas)
    gf = load_gfile_helper(groups_file)
    gf['full_path'] = match_prefix_to_bin(gf['fname'], betas, '.beta')
    return gf


class MarkerFinder:
    def __init__(self, args):
        self.args = args
        self.df = pd.DataFrame()
        self.blocks = pd.DataFrame()
        self.nr_blocks = 0
        self.orig_nr_blocks = 0
        self.keepinds = None
        self.groups = None
        self.verbose = args.verbose
        self.tg_names = []
        self.bg_names = []
        self.first_chunk = {}
        self.hyper, self.hypo = set_hypo_hyper(args.hyper, args.hypo)
        self.validate_args()
        # validate output dir:
        if not op.isdir(args.out_dir):
            os.mkdir(args.out_dir)

        # load groups
        self.gf = load_group_file(args.groups_file, args.betas)
        self.groups = sorted(self.gf['group'].unique())

        # validate target is in groups file
        target = self.args.target
        if target and target not in self.gf['group'].values:
            eprint(f'target {target} not in groups file {self.args.groups_file}')
            eprint('Possible targets:', sorted(self.gf['group'].unique()))
            raise IllegalArgumentError()


    def run(self):

        # load all data
        self.blocks = self.load_blocks_file()
        nr_parts = 20
        chunk_size = int(self.blocks.shape[0] / nr_parts)
        for start in range(0, self.blocks.shape[0], chunk_size):
            self.work_chunck(self.blocks.iloc[start:start + chunk_size])

    def work_chunck(self, blocks_df):
        self.df = self.load_data(blocks_df)

        for group in self.groups:
            if self.args.target and group != self.args.target:
                continue
            eprint(group)
            self.group = group
            self.outpath = op.join(self.args.out_dir, f'Markers.{self.group}.bed')
            if group not in self.first_chunk.keys():
                self.first_chunk[group] = True
            tf = self.find_markers_group()
            self.first_chunk[group] = False
            self.dump_results(tf)

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

    def load_blocks_file(self):
        names = ['chr', 'start', 'end', 'startCpG', 'endCpG']
        cols = range(len(names))
        nrows = DEBUG_NR if self.args.debug else None
        blocks_path = self.args.blocks_path
        df = pd.read_csv(self.args.blocks_path, sep='\t',
                header=None, names=names, nrows=nrows, usecols=cols)
        orig_nr_blocks = df.shape[0]
        df['lenCpG'] = df['endCpG'] - df['startCpG']
        df = df[df['lenCpG'] >= self.args.min_cpg]
        df = df[df['lenCpG'] <= self.args.max_cpg]

        df['len'] = df['end'] - df['start']
        df = df[df['len'] >= self.args.min_bp]
        df = df[df['len'] <= self.args.max_bp]
        df.reset_index(drop=True, inplace=True)

        nr_blocks = df.shape[0]
        if self.verbose:
            eprint(f'loaded {orig_nr_blocks:,}')
            if nr_blocks != orig_nr_blocks:
                eprint(f'droppd to {nr_blocks:,} ')
        return df

    def load_data(self, blocks_df):
        if self.verbose:
            eprint(f'loading {df.shape[0]} blocks...')
        return get_table(blocks_df.copy(), self.gf, self.args.min_cov,
                          self.args.threads, self.verbose, False)

    def set_bg_tg_names(self):
        tg_names = self.gf[self.gf['group'] == self.group]['fname'].values
        bg_names = drop_dup_keep_order(self.gf[self.gf['group'] != self.group]['fname'])

        # remove from bg_names samples shared with tg_names
        bg_names = [s for s in bg_names if s not in tg_names]
        assert(len(bg_names) + len(tg_names) == len(set(self.gf['fname'])))
        return tg_names, bg_names

    def find_markers_group(self):
        if self.df.empty:
            return self.df

        self.tg_names, self.bg_names = self.set_bg_tg_names()

        # look for 'U' markers:
        if not self.hypo:
            tfM = pd.DataFrame()
        else:
            tfU = self.df.copy()
            df_target = tfU[self.tg_names].fillna(.5)
            df_bg = tfU[self.bg_names].fillna(.5)
            tfU[f'tg_low'] = np.quantile(df_target,
                                        self.args.tg_quant,
                                        interpolation='higher',
                                        axis=1)
            tfU[f'bg_high'] = np.quantile(df_bg,
                                         1 - self.args.bg_quant,
                                         interpolation='lower',
                                         axis=1)
            tfU = tfU[(tfU['tg_low'] - tfU['bg_high'] >= self.args.delta)]
            tfU['delta'] = ((tfU['tg_low'] - tfU['bg_high']) * 100).astype(int)
            tfU['direction'] = 'U'

        # look for 'M' markers:
        if not self.hyper:
            tfM = pd.DataFrame()
        else:
            tfM = self.df.copy()
            df_target = tfM[self.tg_names].fillna(.5)
            df_bg = tfM[self.bg_names].fillna(.5)
            tfM[f'tg_high'] = np.quantile(df_target,
                                        1 - self.args.tg_quant,
                                        interpolation='lower',
                                        axis=1)
            tfM[f'bg_low'] = np.quantile(df_bg,
                                         self.args.bg_quant,
                                         interpolation='higher',
                                         axis=1)
            tfM = tfM[(tfM['bg_low'] - tfM['tg_high'] >= self.args.delta)]
            tfM['delta'] = ((tfM['bg_low'] - tfM['tg_high']) * 100).astype(int)
            tfM['direction'] = 'M'
        tf = pd.concat([tfU, tfM]).reset_index(drop=True)
        return tf


    #############################
    #                           #
    #       Dump methods        #
    #                           #
    #############################

    def dump_results(self, tf):
        # eprint(f'Number of markers found: {tf.shape[0]:,}')
        # if self.args.top:
            # tf = tf.head(self.args.top)
        # tf = tf.sort_values(by='startCpG', ascending=True)
        tf['target'] = self.group
        tf['lenCpG'] = tf['lenCpG'].astype(str) + 'CpGs'
        tf['bp'] = (tf['end'] - tf['start']).astype(str) + 'bp'
        tf['region'] = tf['chr'] + ':' + tf['start'].astype(str) + '-' + tf['end'].astype(str)
        tf = tf.round(2)
        cols_to_dump = ['chr', 'start', 'end', 'startCpG', 'endCpG',
                        'target', 'region', 'lenCpG', 'bp', 'tg_low',
                        'bg_high', 'delta', 'direction']
        self.dump_single_df(tf[cols_to_dump], self.outpath)

    def dump_single_df(self, df, outpath):
        """ Dump the results, with comments header """

        eprint(f'dumping to {outpath}')
        df_mode = 'w' if self.first_chunk[self.group] else 'a'
        # write header:
        if self.args.header and self.first_chunk[self.group]:
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
    parser.add_argument('--target', help='find markers only for this group')
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
    # parser.add_argument('--top', type=int,
                        # help='Output only the top TOP markers, under the constraints. [All]')
    parser.add_argument('--header', action='store_true', help='add header to output files')
    parser.add_argument('--debug', '-d', action='store_true', help='Debug mode. Only consider '
                                                                   'first {} blocks'.format(DEBUG_NR))
    parser.add_argument('--tg_quant', type=float, default=.25, help='quantile of target samples to ignore')
    parser.add_argument('--bg_quant', type=float, default=.025, help='quantile of background samples to ignore')
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
