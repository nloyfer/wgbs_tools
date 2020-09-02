#!/usr/bin/python3 -u

import os.path as op
import numpy as np
import pandas as pd
import sys
from utils_wgbs import IllegalArgumentError, eprint, segment_tool, add_GR_args, GenomeRefPaths, \
                       load_dict_section, validate_files_list , validate_single_file, \
                       segment_tool, load_dists, add_multi_thread_args
from genomic_region import GenomicRegion
from multiprocessing import Pool
import argparse
import time
from datetime import timedelta
import subprocess
import multiprocessing

np.set_printoptions(linewidth=200, precision=2)

CHUNK_MAX_SIZE = 60000


def dump_table(df, path):
    df.to_csv(path, sep='\t', header=None, index=None)
    # eprint('dumped to file:', path)


def segment_process(betas, skip, nsites, pcount, max_cpg, max_bp, dist_dict):
    assert nsites, 'trying to segment an empty interval'.format(skip, nsites)
    if nsites == 1:
        return np.array([skip, skip + 1])
    try:
        start_time = time.time()
        beta_files = ' '.join(betas)
        cmd = '{} {} '.format(segment_tool, beta_files)
        cmd += '-s {} -n {} -max_cpg {} '.format(skip, nsites, max_cpg)
        cmd += ' -ps {} -max_bp {} -rd {}'.format(pcount, max_bp, dist_dict)
        brd_str = subprocess.check_output(cmd, shell=True).decode().split()
        # eprint('thread ({}, {}), time: {}'.format(skip, nsites, timedelta(seconds=time.time() - start_time)))
        return np.array(list(map(int, brd_str))) + skip + 1

    except Exception as e:
        eprint('Failed in s={}, n={}'.format(skip, nsites))
        raise e


class SegmentByChunks:
    def __init__(self, args, betas):
        self.gr = GenomicRegion(args)
        self.chunk_size = args.chunk_size
        self.betas = betas
        self.pcount = args.pcount
        self.revdict = self.gr.genome.revdict_path
        self.max_bp = args.max_bp
        self.max_cpg = min(args.max_cpg, args.max_bp // 2)
        assert (self.max_cpg > 1)
        self.args = args

    def break_to_chunks(self):
        """
        Break range of sites to chunks of size 'step'. Examples:

        skip=500, nsites=1500, step=500
        output: sls: [500, 1000]
                nls: [500, 500]

        skip=500, nsites=1550, step=500
        output: sls: [500, 1000]
                nls: [500, 550]
        """
        if self.gr.is_whole():
            nsites = op.getsize(self.betas[0]) // 2
            skip = 0
        else:
            nsites = self.gr.nr_sites
            skip = self.gr.sites[0] - 1

        step = self.chunk_size

        if nsites < step:
            return np.array([skip]), np.array([nsites])

        sls = np.arange(skip, skip + nsites, step)
        nls = np.ones_like(sls) * step

        if len(nls) * step > nsites:  # there's a leftover
            # update the last step:
            nls[-1] = nsites % step

            # in case last step is small, merge is with the one before last
            if nls[-1] < step // 5:
                sls = sls[:-1]
                nls = np.concatenate([nls[:-2], [np.sum(nls[-2:])]])
        return sls, nls

    def run(self):
        skip_list, nsites_list = self.break_to_chunks()
        assert np.all(nsites_list > 0) and len(skip_list) == len(nsites_list) > 0
        p = Pool(self.args.threads)
        params = [(self.betas, si, ni, self.pcount, self.max_cpg, self.max_bp, self.revdict)
                  for si, ni in zip(skip_list, nsites_list)]
        arr = p.starmap(segment_process, params)
        p.close()
        p.join()

        # eprint('merging chunks...')
        df = self.merge_df_list(arr)
        self.dump_result(df)

    def merge_df_list(self, dflist):

        j = 0
        while len(dflist) > 1:
            p = Pool(self.args.threads)
            params = [(dflist[i - 1], dflist[i], self.betas, self.pcount,
                       self.max_cpg, self.max_bp, self.revdict, (j, i // 2 )) for i in range(1, len(dflist), 2)]
            arr = p.starmap(stitch_2_dfs, params)
            p.close()
            p.join()

            last_df = [dflist[-1]] if len(dflist) % 2 else []
            dflist = arr + last_df
            j += 1
        return dflist[0]

    def dump_result(self, df):
        if df is None:
            eprint('Empty blocks array')
            return

        df = np2pd(df)
        df = insert_genomic_loci(df, self.gr)
        dump_table(df, self.args.out_path)


def np2pd(arr):
    return pd.DataFrame({'startCpG': arr[:-1], 'endCpG': arr[1:]})


def stitch_2_dfs(b1, b2, betas, pcount, max_cpg, max_bp, revdict, debind):

    patch1_size = 50
    path2_size = 50
    n1 = b1[-1] - b1[0]
    n2 = b2[-1] - b2[0]
    while patch1_size <= n1 and path2_size <= n2:
        # calculate blocks for patch:
        skip = b1[-1] - patch1_size - 1
        nsites = patch1_size + path2_size
        patch = segment_process(betas, skip, nsites, pcount, max_cpg, max_bp, revdict)

        # find the overlaps
        if is_overlap(b1, patch) and is_overlap(patch, b2):
            # successful stitch with patches 
            return merge2(merge2(b1, patch), b2)
        else:
            # failed stitch - increase patch sizes
            if not is_overlap(b1, patch):
                patch1_size = increase_patch(patch1_size, n1)
            if not is_overlap(patch, b2):
                path2_size = increase_patch(path2_size, n2)

    # Failed: could not stich the two chuncks
    eprint('ERROR: no overlaps at all!!', is_overlap(b1, patch), is_overlap(patch, b2))
    raise IllegalArgumentError('Stitching Failed! Try running with bigger chunk size')


def is_overlap(b1, b2):
    return np.sum(find_dups(b1, b2))


def find_dups(b1, b2):
    return pd.Series(np.concatenate([b1, b2])).duplicated(keep='last')


def merge2(b1, b2):
    nr_from_df1 = np.argmax(np.array(find_dups(b1, b2)))
    skip_from_df2 = np.searchsorted(b2, b1[nr_from_df1])
    return np.concatenate([b1[:nr_from_df1 + 1], b2[skip_from_df2 + 1:]]).copy()


def increase_patch(pre_size, maxval):
    if pre_size == maxval:
        return maxval + 1
    return int(min(pre_size * 2, maxval))


def insert_genomic_loci(df, gr):
    region = None if gr.is_whole() else gr.region_str
    # laod reference dict:
    dict_df = load_dict_section(region, gr.genome_name)

    # merge startCpG
    df = df.merge(dict_df, left_on=['startCpG'], right_on=['idx'])

    # merge endCpG
    dict_df = dict_df.rename(index=str, columns={'idx': 'endCpG', 'start': 'end'})
    df['endCpG'] = df['endCpG'] - 1
    df = df.merge(dict_df[['endCpG', 'end']], on=['endCpG'])
    df['end'] = df['end'] + 1
    df['endCpG'] = df['endCpG'] + 1

    return df[['chr', 'start', 'end', 'startCpG', 'endCpG']]


def parse_args():
    parser = argparse.ArgumentParser()
    add_GR_args(parser)
    betas_or_file = parser.add_mutually_exclusive_group(required=True)
    betas_or_file.add_argument('--betas', nargs='+')
    betas_or_file.add_argument('--beta_file', '-F')
    parser.add_argument('-c', '--chunk_size', type=int, default=CHUNK_MAX_SIZE,
                        help='Chunk size. Default {} sites'.format(CHUNK_MAX_SIZE))
    parser.add_argument('-p', '--pcount', type=float, default=15,
                        help='Pseudo counts of C\'s and T\'s in each block. Default 15')
    parser.add_argument('--max_cpg', type=int, default=1000,
                        help='Maximal allowed blocks size (in #sites). Default is 1000')
    parser.add_argument('--max_bp', type=int, default=2000,
                        help='Maximal allowed blocks size (in bp). Default is 2000')
    parser.add_argument('-o', '--out_path', default=sys.stdout,
                        help='output path [stdout]')
    add_multi_thread_args(parser)
    return parser.parse_args()


def parse_betas_input(args):
    if args.betas:
        betas = args.betas
    elif args.beta_file:
        validate_single_file(args.beta_file)
        with open(args.beta_file, 'r') as f:
            betas = [b.strip() for b in f.readlines() if b.strip() and not b.startswith('#')]
    validate_files_list(betas)
    return betas


def main():
    args = parse_args()
    betas = parse_betas_input(args)
    SegmentByChunks(args, betas).run()


if __name__ == '__main__':
    main()
