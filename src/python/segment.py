#!/usr/bin/python3 -u

import tempfile
import os
import os.path as op
import sys
from multiprocessing import Pool
import argparse
import subprocess
import numpy as np
import pandas as pd
from utils_wgbs import IllegalArgumentError, eprint, segment_tool, add_GR_args, \
                       validate_file_list, validate_single_file, \
                       add_multi_thread_args, GenomeRefPaths, validate_local_exe, \
                       beta_sanity_check
from convert import add_bed_to_cpgs
from genomic_region import GenomicRegion, index2chrom
from beta_to_blocks import load_blocks_file


DEF_CHUNK = 60000


# TODO: remove it?
def is_block_file_nice(df):

    # no duplicated blocks
    if df.shape[0] != df.drop_duplicates().shape[0]:
        msg = 'Some blocks are duplicated'
        return False, msg

    # no overlaps between blocks
    sdf = df.sort_values(by='startCpG')
    if not (sdf['startCpG'][1:].values - sdf['endCpG'][:sdf.shape[0] - 1].values  >= 0).all():
        msg = 'Some blocks overlap'
        return False, msg

    return True, ''


def segment_process(params):
    sites = params['sites']
    start, end = sites
    assert end - start > 0, f'trying to segment an empty interval {sites}'
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
        self.betas = betas
        max_cpg = min(args.max_cpg, args.max_bp // 2)
        assert max_cpg > 1
        self.genome = GenomeRefPaths(args.genome)
        self.param_dict = {'betas': betas,
                          'pcount': args.pcount,
                          'max_cpg': max_cpg,
                          'max_bp': args.max_bp,
                          'revdict': self.genome.revdict_path,
                          'genome': self.genome
                          }
        self.args = args
        self.validate_genome()

    def validate_genome(self):
        for beta in self.betas:
            if not beta_sanity_check(beta, self.genome):
                msg = f'[wt segment] ERROR: current genome reference ({self.genome.genome}) does not match the input beta file ({beta}).'
                raise IllegalArgumentError(msg)

    def break_to_chunks(self):
        """ Break range of sites to chunks of size 'step',
            while keeping chromosomes separated """
        # print a warning in case chunk size is too small
        step = self.args.chunk_size
        if step < self.args.max_cpg:
            msg = '[wt segment] WARNING: chunk_size is small compared to max_cpg and/or max_bp.\n' \
                  '                      It may cause wt segment to fail. It\'s best setting\n' \
                  '                      chunk_size > min{max_cpg, max_bp/2}'
            eprint(msg)

        if self.args.bed_file:
            df = load_blocks_file(self.args.bed_file)[['startCpG', 'endCpG']].dropna()
            # make sure bed file has no overlaps or duplicated regions
            is_nice, msg = is_block_file_nice(df)
            if not is_nice:
                msg = '[wt segment] ERROR: invalid bed file.\n' \
                      f'                    {msg}\n' \
                      f'                    Try: sort -k1,1 -k2,2n {self.args.bed_file} | ' \
                      'bedtools merge -i - | wgbstools convert --drop_empty -p -L -'
                eprint(msg)
                raise IllegalArgumentError('Invalid bed file')
            if df.shape[0] > 2*1e4:
                msg = '[wt segment] WARNING: bed file contains many regions.\n' \
                      '                      Segmentation will take a long time.\n' \
                      '                      Consider running w/o -L flag and intersect the results\n'
                eprint(msg)

        else:   # No bed file provided
            gr = GenomicRegion(self.args)
            # whole genome - make a dummy "bed file" of the full chromosomes
            if gr.is_whole():
                cf = self.genome.get_chrom_cpg_size_table()
                cf['endCpG'] = np.cumsum(cf['size']) + 1
                cf['startCpG'] = cf['endCpG'] - cf['size']
                df = cf[['startCpG', 'endCpG']]
            # one region
            else:
                df = pd.DataFrame(columns=['startCpG', 'endCpG'], data=[gr.sites])

        # build a DataFrame of chunks, with a "tag"/label field,
        # so we know which chunks to merge later on.
        tags = []
        starts = []
        ends = []
        for _, row in df.iterrows():
            start, end = row
            bords = list(range(start, end, step)) + [end]
            tags += [f'{start}-{end}'] * (len(bords) -1)
            starts += bords[:-1]
            ends += bords[1:]
        return tags, starts, ends

    def run(self):
        # break input region/s to small chunks
        tags, starts, ends = self.break_to_chunks()
        # segment each chunk separately in a single thread
        p = Pool(self.args.threads)
        params = [(dict(self.param_dict, **{'sites': (s, e)}),) for s, e in zip(starts, ends)]
        arr = p.starmap(segment_process, params)
        p.close()
        p.join()

        # merge chunks from the same "tag" group
        # (i.e. the same chromosome, or the same region of the provided bed file)
        df = pd.DataFrame()
        for tag in set(tags):
            carr = [arr[i] for i in range(len(arr)) if tags[i] == tag]
            merged = self.merge_df_list(carr)
            df = pd.concat([df, pd.DataFrame({'startCpG': merged[:-1], 'endCpG': merged[1:]})])
        self.dump_result(df.reset_index(drop=True))

    def merge_df_list(self, dflist):
        # Given a set of chunks to merge, recursively pairwise stich them.

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
        if df.empty:
            eprint('Empty blocks array')
            return

        # sort by startCpG and filter by CpGs
        nr_blocks = df.shape[0]
        df.sort_values(by=['startCpG'], inplace=True)
        df = df[df.endCpG - df.startCpG > self.args.min_cpg - 1].reset_index(drop=True)

        # verbose
        nr_blocks_filt = df.shape[0]
        nr_dropped = nr_blocks - nr_blocks_filt
        eprint(f'[wt segment] found {nr_blocks_filt:,} blocks\n' \
               f'             (dropped {nr_dropped:,} short blocks)')

        # add genomic loci and dump/print
        temp_path = next(tempfile._get_candidate_names())
        try:
            df.to_csv(temp_path, sep='\t', header=None, index=None)
            add_bed_to_cpgs(temp_path, self.genome.genome, self.args.out_path)
        finally:
            if op.isfile(temp_path):
                os.remove(temp_path)


#############################################################
#                                                           #
#           Chunk stiching logic                            #
#                                                           #
#############################################################

def stitch_2_dfs(b1, b2, params):

    # if b2 is not the direct extension of b1, we have a problem
    if b1[-1] != b2[0]:
        msg = '[wt segment] Patch stitching Failed! ' \
              '             patches are not supposed to be merged'
        raise IllegalArgumentError(msg)

    n1 = b1[-1] - b1[0]
    n2 = b2[-1] - b2[0]
    patch1_size = min(50, n1)
    patch2_size = min(50, n2)
    patch = np.array([], dtype=int)
    while patch1_size <= n1 and patch2_size <= n2:
        # calculate blocks for patch:
        start = b1[-1] - patch1_size #- 1
        end = b1[-1] + patch2_size
        cparams = dict(params, **{'sites': (start, end)})
        patch = segment_process(cparams)

        # find the overlaps
        if is_2_overlap(b1, patch) and is_2_overlap(patch, b2):
            # successful stitch with patches
            return merge2(merge2(b1, patch), b2)
        else:
            # failed stitch - increase patch sizes
            if not is_2_overlap(b1, patch):
                patch1_size = increase_patch(patch1_size, n1)
            if not is_2_overlap(patch, b2):
                patch2_size = increase_patch(patch2_size, n2)

    # Failed: could not stich the two chuncks
    msg = '[wt segment] Patch stitching Failed! ' \
          '             Try increasing chunk size (--chunk_size flag)'
    raise IllegalArgumentError(msg)


def is_2_overlap(b1, b2):
    return np.sum(find_dups(b1, b2))


def find_dups(b1, b2):
    return pd.Series(np.concatenate([b1, b2])).duplicated(keep=False).values


def merge2(b1, b2):
    nr_from_df1 = np.argmax(find_dups(b1, b2))
    skip_from_df2 = np.searchsorted(b2, b1[nr_from_df1])
    return np.concatenate([b1[:nr_from_df1 + 1], b2[skip_from_df2 + 1:]]).copy()


def increase_patch(pre_size, maxval):
    if pre_size == maxval:
        return maxval + 1  # too large, so the while loop will break
    return int(min(pre_size * 2, maxval))


#############################################################
#                                                           #
#                       Main                                #
#                                                           #
#############################################################

def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    add_GR_args(parser, bed_file=True)
    betas_or_file = parser.add_mutually_exclusive_group(required=True)
    betas_or_file.add_argument('--betas', nargs='+')
    betas_or_file.add_argument('--beta_file', '-F')
    parser.add_argument('-c', '--chunk_size', type=int, default=DEF_CHUNK,
                        help=f'Chunk size. Default {DEF_CHUNK} sites')
    parser.add_argument('-p', '--pcount', type=float, default=15,
                        help='Pseudo counts of C\'s and T\'s in each block. Default 15')
    parser.add_argument('--min_cpg', type=int, default=1,
                        help='Minimal block size (in #sites) to output. Shorter blocks will simply be ' \
                             'ommited from output (equivalent to set min_cpg to 1 and then filter output by ' \
                             'length). Default is 1')
    parser.add_argument('--max_cpg', type=int, default=1000,
                        help='Maximal allowed block size (in #sites). Default is 1000')
    parser.add_argument('--max_bp', type=int, default=2000,
                        help='Maximal allowed block size (in bp). Default is 2000')
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
    validate_local_exe(segment_tool)
    betas = parse_betas_input(args)
    SegmentByChunks(args, betas).run()


if __name__ == '__main__':
    main()
