import os.path as op
import numpy as np
import pandas as pd
import sys
from utils_wgbs import IllegalArgumentError, eprint, pat_segment_tool, \
    load_dict_section, segment_tool, GenomeRefPaths
from genomic_region import GenomicRegion
from multiprocessing import Pool
from pat_processor_client import PatProcesserOutput
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


def get_pat_processor_command(skip, nr_sites):
    start = skip
    end = skip + nr_sites
    sites_str = "{}-{}".format(start, end)
    try:
        GenomicRegion()
    except Exception:
        x=0


def get_segmentation_command(x, data_files, skip, nsites, max_block_size, pcount):
    return {
        PatProcesserOutput.BETA: '{} {} -s {} -n {} -m {} -ps {} -st {}'.format(segment_tool, data_files, skip,
                                                      nsites, max_block_size, pcount, "beta"),
        PatProcesserOutput.RHO: '{} {} -s {} -n {} -m {} -ps {} -st {}'.format(segment_tool, data_files, skip,
                                                      nsites, max_block_size, pcount, "rho"),
        PatProcesserOutput.LBETA: '{} {} -s {} -n {} -m {} -ps {} -st {}'.format(segment_tool, data_files, skip,
                                                                               nsites, max_block_size, pcount, "lbeta"),
        PatProcesserOutput.PAT: '{} {} -s {} -n {} -m {} -ps {}'.format(pat_segment_tool, data_files, skip,
                                                                               nsites, max_block_size, pcount)
    }[x]


def segment_process(betas, skip, nsites, pcount, max_block_size, pat_processor_output_type):
    assert nsites, 'trying to segment an empty interval'.format(skip, nsites)
    if nsites == 1:
        return np.array([skip, skip + 1])
    try:
        start_time = time.time()
        beta_files = ' '.join(betas)
        cmd = get_segmentation_command(pat_processor_output_type, beta_files, skip, nsites, max_block_size, pcount)
        brd_str = subprocess.check_output(cmd, shell=True).decode().split()
        # eprint('thread ({}, {}), time: {}'.format(skip, nsites, timedelta(seconds=time.time() - start_time)))
        return np.array(list(map(int, brd_str))) + skip + 1

    except Exception as e:
        eprint('Failed in s={}, n={}'.format(skip, nsites))
        raise e


class SegmentByChunks:
    def __init__(self, args, in_files, pat_processor_output_type):
        self.gr = GenomicRegion(args)
        self.chunk_size = args.chunk_size
        self.in_files = in_files
        if args.max_bp:
            self.max_block_size = args.max_bp
        if args.pcount:
            self.pcount = args.pcount
        self.args = args
        self.pat_processor_output_type = pat_processor_output_type

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
            nsites = op.getsize(self.in_files[0]) // 2
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

    def break_to_chunks_by_chrome(self):
        """
        Break range of sites to chunks of size 'step'. Examples:

        skip=500, nsites=1500, step=500
        output: sls: [500, 1000]
                nls: [500, 500]

        skip=500, nsites=1550, step=500
        output: sls: [500, 1000]
                nls: [500, 550]
        """
        gr = GenomicRegion(region="chr1")
        chrom_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                      "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                      "chr22", "chrX", "chrY", "chrM"]
        cpg_file_path = GenomeRefPaths().dict_path
        chrom_dict = {}
        for chrom in chrom_list:
            tabix_head_cmd = "tabix {} {} | head -1".format(cpg_file_path, chrom)
            tabix_tail_cmd = "tabix {} {} | tail -1".format(cpg_file_path, chrom)
            head_str = subprocess.check_output(tabix_head_cmd, shell=True).decode().split()
            tail_str = subprocess.check_output(tabix_tail_cmd, shell=True).decode().split()
            start_ind = int(head_str[-1])
            end_ind = int(tail_str[-1])
            x = 0
        if self.gr.is_whole():
            nsites = op.getsize(self.in_files[0]) // 2
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
        params = [(self.in_files, si, ni, self.pcount, self.max_block_size, self.pat_processor_output_type)
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
            params = [(dflist[i - 1], dflist[i], self.in_files, self.pcount,
                       self.max_block_size, (j, i // 2), self.pat_processor_output_type) for i in range(1, len(dflist), 2)]
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

    def dump_pat_blocks(self, blocks_list):
        if blocks_list is None:
            eprint('Empty blocks array')
            return

        flat_list = [(item.split(" ")[0], item.split(" ")[1]) for item in blocks_list]
        # chr_arr = np.array([chrome for chrome, _, _ in flat_list])
        start_arr = np.array([int(start_ind) for start_ind, _ in flat_list])
        end_arr = np.array([int(end_ind) for _, end_ind in flat_list])
        df = pd.DataFrame({'startCpG': start_arr, 'endCpG': end_arr})
        df = insert_genomic_loci(df, self.gr)

        dump_table(df, self.args.out_path)

    def dump_cpg_counts(self, blocks_list):
        if blocks_list is None:
            eprint('Empty blocks array')
            return

        flat_list = [item.split("\t") for item in blocks_list]
        # chr_arr = np.array([chrome for chrome, _, _ in flat_list])
        end_arr = np.array([int(el[0]) for el in flat_list])
        start_arr = end_arr
        u_arr = np.array([int(el[1]) for el in flat_list])
        x_arr = np.array([int(el[2]) for el in flat_list])
        m_arr = np.array([int(el[3]) for el in flat_list])
        end_arr = end_arr + 1
        df = pd.DataFrame({'startCpG': start_arr, 'endCpG': end_arr})
        df = insert_genomic_loci(df, self.gr)
        df['u'] = u_arr
        df['x'] = x_arr
        df['m'] = m_arr
        dump_table(df, self.args.out_path)


def np2pd(arr):
    return pd.DataFrame({'startCpG': arr[:-1], 'endCpG': arr[1:]})


def stitch_2_dfs(b1, b2, betas, pcount, max_block_size, debind, pat_processor_output_type):

    patch1_size = 50
    path2_size = 50
    n1 = b1[-1] - b1[0]
    n2 = b2[-1] - b2[0]
    while patch1_size <= n1 and path2_size <= n2:
        # calculate blocks for patch:
        skip = b1[-1] - patch1_size - 1
        nsites = patch1_size + path2_size
        patch = segment_process(betas, skip, nsites, pcount, max_block_size, pat_processor_output_type)

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