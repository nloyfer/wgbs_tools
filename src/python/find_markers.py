#!/usr/bin/python3 -u

import os
import os.path as op
import subprocess
import sys
import numpy as np
import pandas as pd

from beta_to_blocks import load_blocks_file
from betas_to_table import get_table, groups_load_wrap
from dmb import load_gfile_helper, match_prefix_to_bin
from utils_wgbs import validate_single_file, add_multi_thread_args, \
        validate_file_list, eprint, drop_dup_keep_order, IllegalArgumentError
from fm_load_params import MFParams, parse_args



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

        # validate output dir:
        if not op.isdir(args.out_dir):
            os.mkdir(args.out_dir)

        # load groups
        self.gf = load_group_file(args.groups_file, args.betas)
        groups = sorted(self.gf['group'].unique())
        self.targets = get_validate_targets(self.args.targets, groups)
        self.res = {t: pd.DataFrame() for t in self.targets}
        self.inds_dict = set_bg_tg_names(self.gf, self.targets)

    def run(self):

        # first dump the param file
        self.dump_params()

        # load all blocks
        blocks = self.load_blocks()
        # load data in chunks
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
            tf = self.find_group_markers()
            self.res[group] = pd.concat([self.res[group], tf])

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

    def find_group_markers(self):
        if self.df.empty:
            return self.df

        # load context
        self.tg_names, self.bg_names = self.inds_dict[self.group]
        tf = self.df.copy()

        # filter blocks by coverage:
        df_tg = tf[self.tg_names]
        df_bg = tf[self.bg_names]
        keep_tg = (df_tg.notna().sum(axis=1) / (df_tg.shape[1])) >= (1 - self.args.na_rate_tg)
        keep_bg = (df_bg.notna().sum(axis=1) / (df_bg.shape[1])) >= (1 - self.args.na_rate_bg)
        tf = tf.loc[keep_tg & keep_bg, :].reset_index(drop=True)

        # find markers:
        tfM = self.find_M_markers(tf)
        tfU = self.find_U_markers(tf)
        tf = pd.concat([tfU, tfM]).reset_index(drop=True)
        tf['delta'] = ((tf['tg'] - tf['bg']).abs() * 100).astype(int)
        return tf

    def find_M_markers(self, tf):
        # look for 'M' markers
        if self.args.only_hypo:
            return pd.DataFrame()

        # filter by mean thresholds
        df_bg = tf[self.bg_names]
        df_tg = tf[self.tg_names]
        keep_umt = np.nanmean(df_bg.values, axis=1) <= self.args.unmeth_mean_thresh
        keep_mmt = np.nanmean(df_tg.values, axis=1) >= self.args.meth_mean_thresh
        keep = keep_umt & keep_mmt

        tfM = tf.loc[keep, :].copy()
        tfM['direction'] = 'M'
        df_tg = tf.loc[keep, self.tg_names].fillna(.5)
        df_bg = tf.loc[keep, self.bg_names].fillna(.5)

        # filter by quantile thresholds
        tfM['tg'] = np.quantile(df_tg,
                                self.args.tg_quant,
                                interpolation='higher',
                                axis=1)
        tfM['bg'] = np.quantile(df_bg,
                                1 - self.args.bg_quant,
                                interpolation='lower',
                                axis=1)
        # delta
        tfM = tfM[(tfM['tg'] - tfM['bg'] >= self.args.delta)].reset_index(drop=True)
        if tfM.empty:
            return tfM
        # high & low
        keep_ut = tfM['bg'] <= self.args.unmeth_thresh
        keep_mt = tfM['tg'] >= self.args.meth_thresh
        keep = keep_ut & keep_mt
        tfM = tfM.loc[keep, :]

        return tfM

    def find_U_markers(self, tf):
        # look for 'U' markers
        if self.args.only_hyper:
            return pd.DataFrame()

        df_bg = tf[self.bg_names]
        df_tg = tf[self.tg_names]
        keep_umt = np.nanmean(df_tg.values, axis=1) <= self.args.unmeth_mean_thresh
        keep_mmt = np.nanmean(df_bg.values, axis=1) >= self.args.meth_mean_thresh
        keep = keep_umt & keep_mmt

        tfU = tf.loc[keep, :].copy()
        tfU['direction'] = 'U'
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
        # delta
        tfU = tfU[(tfU['bg'] - tfU['tg'] >= self.args.delta)].reset_index(drop=True)
        if tfU.empty:
            return tfU
        # high & low
        keep_ut = tfU['tg'] <= self.args.unmeth_thresh
        keep_mt = tfU['bg'] >= self.args.meth_thresh
        keep = keep_ut & keep_mt
        tfU = tfU.loc[keep, :]
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

        # todo: add all parameters to a config file in the same directory. or to header
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

    def dump_params(self):
        """ Dump a parameter file """
        outpath = op.join(self.args.out_dir, 'params.txt')
        with open(outpath, 'w') as f:
            for key in vars(self.args):
                val = getattr(self.args, key)
                if key == 'beta_list_file':
                    val = None
                if key == 'betas':
                    val = ' '.join(val)
                f.write(f'{key}:{val}\n' )
                # f.write(f'#> {sample}\n' )
        eprint(f'dumped parameter file to {outpath}')


def main():
    """
    Find differentially methylated blocks
    """
    params = MFParams(parse_args())
    MarkerFinder(params).run()


if __name__ == '__main__':
    main()
