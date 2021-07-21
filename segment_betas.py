#!/usr/bin/python3 -u

import os.path as op
import numpy as np
import pandas as pd
import sys
from utils_wgbs import IllegalArgumentError, eprint, segment_tool, add_GR_args, \
                       load_dict_section, validate_file_list , validate_single_file, \
                       add_multi_thread_args, GenomeRefPaths
from genomic_region import GenomicRegion, index2chrom
from multiprocessing import Pool
import argparse
import subprocess
import multiprocessing


DEF_CHUNK = 60000

def break_to_chunks_helper(start, end, step):
    res = []
    while start + step < end:
        res.append((start, start + step))
        start = start + step
    res.append((start, end))
    return res


def segment_process(params):
    sites = params['sites']
    start, end = sites
    assert end - start >= 1, f'trying to segment an empty interval {sites}'
    if end - start == 1:
        return np.array([start, end])
    try:
        beta_files = ' '.join(params['betas'])
        cmd = f'{segment_tool} {beta_files} '
        cmd += f'-s {start - 1} -n {end - start} -max_cpg {params["max_cpg"]} '
        cmd += f' -ps {params["pcount"]} -max_bp {params["max_bp"]} '
        chrom = index2chrom(start, params["genome"])
        cmd = f'tabix {params["revdict"]} {chrom}:{start}-{end - 1} | cut -f2 |' + cmd
        brd_str = subprocess.check_output(cmd, shell=True).decode().split()
        return np.array(list(map(int, brd_str))) + start

    except Exception as e:
        eprint(f'Failed in sites {sites}')
        raise e


class SegmentByChunks:
    def __init__(self, args, betas):
        self.gr = GenomicRegion(args)       # TODO: support -L (e.g. for MCC-seq or TWIST array)
        self.betas = betas
        max_cpg = min(args.max_cpg, args.max_bp // 2)
        assert (max_cpg > 1)
        self.genome = GenomeRefPaths(args.genome)
        self.param_dict = {'betas': betas,
                          'pcount': args.pcount,
                          'max_cpg': max_cpg,
                          'max_bp': args.max_bp,
                          'revdict': self.genome.revdict_path,
                          'genome': self.genome
                          }
        self.args = args
        if args.chunk_size < max_cpg:
            msg = '[wt segment] WARNING: chunk_size is small compared to max_cpg and/or max_bp.\n' \
                  '                      It may cause wt segment to fail. It\'s best setting\n' \
                  '                      chunk_size > min{max_cpg, max_bp/2}'
            eprint(msg)

    def break_to_chunks(self):
        """ Break range of sites to chunks of size 'step', while keeping chromosomes separated """
        step = self.args.chunk_size
        if not self.gr.is_whole():
            start, end = self.gr.sites
            return break_to_chunks_helper(start, end, step)
        else:
            res = []
            cf = self.genome.get_chrom_cpg_size_table()
            cf['borders'] = np.cumsum(cf['size'])
            for _, row in cf.iterrows():
                chrom, size, border = row
                r = break_to_chunks_helper(border - size + 1, border + 1, step)
                res += r
        return res

    def run(self):
        chunks = self.break_to_chunks()
        # print(np.array(chunks))
        p = Pool(self.args.threads)
        params = [(dict(self.param_dict, **{'sites': s}),) for s in chunks]
        arr = p.starmap(segment_process, params)
        p.close()
        p.join()

        df = self.merge_df_list(arr)
        self.dump_result(df)

    def merge_df_list(self, dflist):

        while len(dflist) > 1:
            p = Pool(self.args.threads)
            params = [(dflist[i - 1], dflist[i], self.param_dict) for i in range(1, len(dflist), 2)]
            arr = p.starmap(stitch_2_dfs, params)
            p.close()
            p.join()

            last_df = [dflist[-1]] if len(dflist) % 2 else []
            dflist = arr + last_df
        return dflist[0]

    def dump_result(self, df):
        if df is None:
            eprint('Empty blocks array')
            return

        df = pd.DataFrame({'startCpG': df[:-1], 'endCpG': df[1:]})
        eprint(f'[wt segment] found {df.shape[0]} blocks')
        df = insert_genomic_loci(df, self.gr)
        df.to_csv(self.args.out_path, sep='\t', header=None, index=None)


def stitch_2_dfs(b1, b2, params):

    # if b1 and b2 are the edges of different chromosomes, simply concatenate them
    if index2chrom(b1[-1] - 1, params['genome']) != index2chrom(b2[0], params['genome']):
        return np.concatenate([b1[:-1], b2])

    n1 = b1[-1] - b1[0]
    n2 = b2[-1] - b2[0]
    patch1_size = min(50, n1)
    patch2_size = min(50, n2)
    patch = np.array([], dtype=int)
    while patch1_size <= n1 and patch2_size <= n2:
        # calculate blocks for patch:
        start = b1[-1] - patch1_size - 1
        end = b1[-1] + patch2_size
        cparams = dict(params, **{'sites': (start, end)})
        patch = segment_process(cparams)

        # find the overlaps
        if is_overlap(b1, patch) and is_overlap(patch, b2):
            # successful stitch with patches 
            return merge2(merge2(b1, patch), b2)
        else:
            # failed stitch - increase patch sizes
            if not is_overlap(b1, patch):
                patch1_size = increase_patch(patch1_size, n1)
            if not is_overlap(patch, b2):
                patch2_size = increase_patch(patch2_size, n2)

    # Failed: could not stich the two chuncks
    # eprint('ERROR: no overlaps at all!!', is_overlap(b1, patch), is_overlap(patch, b2))
    msg = '[wt segment] Patch stitching Failed! Try increasing chunk size (--chunk_size flag)'
    raise IllegalArgumentError(msg)


def is_overlap(b1, b2):
    return np.sum(find_dups(b1, b2))


def find_dups(b1, b2):
    return pd.Series(np.concatenate([b1, b2])).duplicated(keep=False).values


def merge2(b1, b2):
    nr_from_df1 = np.argmax(find_dups(b1, b2))
    skip_from_df2 = np.searchsorted(b2, b1[nr_from_df1])
    return np.concatenate([b1[:nr_from_df1 + 1], b2[skip_from_df2 + 1:]]).copy()


def increase_patch(pre_size, maxval):
    # TODO: limit patch to not cross chromosomes.
    # Currently the c++ tool will throw exception if it happens.
    if pre_size == maxval:
        return maxval + 1  # too large, so the while loop will break
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
    parser = argparse.ArgumentParser(description=main.__doc__)
    add_GR_args(parser)
    betas_or_file = parser.add_mutually_exclusive_group(required=True)
    betas_or_file.add_argument('--betas', nargs='+')
    betas_or_file.add_argument('--beta_file', '-F')
    parser.add_argument('-c', '--chunk_size', type=int, default=DEF_CHUNK,
                        help=f'Chunk size. Default {DEF_CHUNK} sites')
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
    """
    parse user input to get the list of beta files to segment
    Either args.betas is a list of beta files,
    or args.beta_file is a text file in which each line is a beta file
    return: list of beta files
    """
    if args.betas:
        betas = args.betas
    elif args.beta_file:
        validate_single_file(args.beta_file)
        with open(args.beta_file, 'r') as f:
            betas = [b.strip() for b in f.readlines() if b.strip() and not b.startswith('#')]
        if not betas:
            raise IllegalArgumentError(f'no beta files found in file {args.beta_file}')
    validate_file_list(betas)
    return betas


def main():
    """
    Segment the genome, or a subset region, to homogenously methylated blocks.
    Input: one or more beta files to segment
    Output: blocks file (BED format + startCpG, endCpG columns)
    """
    args = parse_args()
    betas = parse_betas_input(args)
    SegmentByChunks(args, betas).run()


if __name__ == '__main__':
    main()
