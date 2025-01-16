import argparse
import subprocess
import sys
import warnings
import re
from genomic_region import GenomicRegion
from utils_wgbs import add_GR_args, MAX_PAT_LEN, GenomeRefPaths, add_multi_thread_args, \
    IllegalArgumentError, eprint, COORDS_COLS5
from multiprocessing import Pool
import pandas as pd
try:
    from statsmodels.stats.multitest import multipletests
except:
    raise IllegalArgumentError("Please install statsmodels package in order to use this feature. i.e. 'pip install statsmodels'")

import numpy as np
from scipy import stats


c_char = b'C'
t_char = b'T'
# c_char = 'C'
# t_char = 'T'

def read_pat_vis(pat_text, start, end, is_strict, min_len):
    first_ind = 0
    max_ind = 0

    reads_list = []
    ind_list = []
    counts_list = []
    is_first = True
    for line in pat_text.split("\n"):
        # line = line.strip()
        if len(line) > 1:
            tokens = line.split("\t")
            cur_ind = int(tokens[1])
            pat = tokens[-2]
            cur_end = cur_ind + len(pat)
            if cur_end <= start:
                continue
            if is_strict:
                if cur_ind < start:
                    pat = pat[start - cur_ind:]
                    cur_ind = start
                if cur_ind + len(pat) > end:
                    pat = pat[:end - cur_ind]
            if len(pat) < min_len:
                continue

            if is_first:
                first_ind = cur_ind
                is_first = False
            last_ind = cur_ind
            ind_list.append(last_ind)
            reads_list.append(pat)
            counts_list.append(int(tokens[-1]))
            if cur_end > max_ind:
                max_ind = cur_end
    num_cpgs = max_ind - first_ind
    pat_matrix = np.chararray((sum(counts_list), num_cpgs))
    pat_matrix[:] = b'.'
    read_index = 0
    for ((read, read_ind), count) in zip(zip(reads_list, ind_list), counts_list):
        col_ind = read_ind - first_ind
        for _ in range(count):
            pat_matrix[read_index, col_ind:(col_ind + len(read))] = list(read)
            read_index += 1
    return pat_matrix


def em_pat_matrix(pat_matrix, should_print=True):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        c_per_col = 1e-3 + sum(pat_matrix == c_char)
        t_per_col = 1e-3 + sum(pat_matrix == t_char)
        n_per_col = c_per_col + t_per_col

        num_reads = pat_matrix.shape[0]
        num_cpgs = pat_matrix.shape[1]

        p_c = np.zeros([2, num_cpgs])
        p_c[0, :] = 0.9
        p_c[1, :] = 0.1
        p_t = 1 - p_c
        l_p_c = np.log2(p_c)
        l_p_t = np.log2(p_t)
        p_alleles = np.array([0.5, 0.5])
        l_p_alleles = np.log2(p_alleles)
        ll = np.inf * -1

        ll_alleles = np.zeros([2, num_reads])
        ll_delta = 100
        while ll_delta > 0:
            # for row_ind in range(num_reads):
            #     #l_p_alleles[0]*np.ones(num_reads) + np.matmul(l_p_c[0,:], (pat_matrix == c_char).T) + np.matmul(l_p_t[0, :], (pat_matrix == t_char).T)
            #     ll_alleles_c[0, row_ind] = l_p_alleles[0] + sum(l_p_c[0, pat_matrix[row_ind, :] == c_char]) \
            #                              + sum(l_p_t[0, pat_matrix[row_ind, :] == t_char])
            #     ll_alleles_c[1, row_ind] = l_p_alleles[1] + sum(l_p_c[1, pat_matrix[row_ind, :] == c_char]) + \
            #                              sum(l_p_t[1, pat_matrix[row_ind, :] == t_char])
            ll_alleles[0, :] = l_p_alleles[0] * np.ones(num_reads) + \
                               np.matmul(l_p_c[0, :], (pat_matrix == c_char).T) + \
                               np.matmul(l_p_t[0, :], (pat_matrix == t_char).T)
            ll_alleles[1, :] = l_p_alleles[1] * np.ones(num_reads) + \
                               np.matmul(l_p_c[1, :], (pat_matrix == c_char).T) + \
                               np.matmul(l_p_t[1, :], (pat_matrix == t_char).T)
            # read_assignments_o = np.argmax(ll_alleles_c, axis=0)
            read_assignments = np.argmax(ll_alleles, axis=0)
            new_ll = sum(ll_alleles[0, read_assignments == 0]) + sum(ll_alleles[1, read_assignments == 1])
            ll_delta = new_ll - ll
            ll = new_ll

            # p_alleles = np.array([sum(read_assignments == 0), sum(read_assignments == 1)]) / num_reads
            # l_p_alleles = np.log2(p_alleles)
            # for j in range(2):
            #     for i in range(num_cpgs):
            #         p_c[j, i] = 1e-3 + sum(pat_matrix[read_assignments == j, i] == c_char)
            #         p_t[j, i] = 1e-3 + sum(pat_matrix[read_assignments == j, i] == t_char)
            for z in range(2):
                p_c[z, :] = 1e-3 + sum(pat_matrix[read_assignments == z, :] == c_char)
                p_t[z, :] = 1e-3 + sum(pat_matrix[read_assignments == z, :] == t_char)
            totals = p_c + p_t
            l_p_c = np.log2(p_c / totals)
            l_p_t = np.log2(p_t / totals)
        if should_print:
            bpi = 2 ** (ll / sum(n_per_col))
            print(f"LL: {ll} | {num_reads} reads | {int(round(sum(n_per_col)))} observed | BPI: {bpi}")
        return ll


def calc_initial_liklihood(pat_matrix, should_print=True):
    try:
        c_per_col = 1e-3 + sum(pat_matrix == c_char)
        t_per_col = 1e-3 + sum(pat_matrix == t_char)
        n_per_col = c_per_col + t_per_col
        l_p_c = np.log2(c_per_col / n_per_col)
        l_p_t = np.log2(t_per_col / n_per_col)
        ll0 = sum(sum(pat_matrix == c_char) * l_p_c + sum(pat_matrix == t_char) * l_p_t)
        bpi = 2 ** (ll0 / sum(n_per_col))
        if should_print:
            print(f"LL: {ll0} | {pat_matrix.shape[0]} reads | {int(round(sum(n_per_col)))} observed | BPI: {bpi}")
        return ll0
    except:
        return 0


def pull_pat_file(region, file_name):
    cmd = f'tabix {file_name} {region}'
    return subprocess.check_output(cmd, shell=True).decode()


def test_single_region(pat_file, chrom, sites, is_strict, min_len, should_print=True):
    s1, s2 = sites
    new_start = max(1, s1 - MAX_PAT_LEN)
    pat_region = f"{chrom}:{new_start}-{s2-1}"

    # pat_parser = vis.parse_args()
    # pat_args_list = [pat_file, "-r", region, "--no_dense", "--no_color", "--max_reps", "999", "--text"]
    # pat_args = pat_parser.parse_args(pat_args_list)
    #
    #
    # pat_mat = PatVis(pat_args, pat_file).get_block()["table"]
    pat_text = pull_pat_file(pat_region, pat_file)
    pat_mat = read_pat_vis(pat_text, s1, s2, is_strict, min_len)
    if pat_mat.shape[0] == 0:
        return np.float32(1.0)

    ll0 = calc_initial_liklihood(pat_mat, should_print=should_print)
    post_em_ll = em_pat_matrix(pat_mat, should_print=should_print)
    test_stat = 2 * np.log(2) * (post_em_ll - ll0)
    pv = 1 - stats.chi2.cdf(test_stat, pat_mat.shape[1])
    if should_print:
        print(f"pvalue: {pv:,.3e}")
    return pv


def read_blocks_and_test(tabixed_bed_file, cur_region, pat_file, is_strict, min_len, verbose=False):
    tabix_cmd = f"tabix {tabixed_bed_file} {cur_region}"
    cur_blocks_lines = subprocess.check_output(tabix_cmd, shell=True).decode().split("\n")
    p_val_list = []
    for line in cur_blocks_lines:
        if not line.strip():
            continue
        tokens = line.split("\t")
        sites = (int(tokens[3]), int(tokens[4]))
        p_val = test_single_region(pat_file, tokens[0], sites, is_strict, min_len, should_print=False)
        p_val = p_val.astype(np.float32)
        p_val_list.append((line, p_val))
    if verbose:
        eprint(f"[wt bimodal] finished processesing {cur_region}")
    return p_val_list


def choose_blocks_by_fdr_bh(pvals, blocks, alpha=0.05, return_all=False):
    rejected_list, corrected_p_vals, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')
    if not rejected_list[0]:
        return [], []
    # keep only the first "index" values (the significant ones)
    if not return_all:
        index = np.argmax(1 - rejected_list)
        return blocks[:index], corrected_p_vals[:index]
    else:
        return blocks, corrected_p_vals


def test_multiple_regions(tabixed_bed_file, pat_file, num_threads, out_file, is_strict, min_len, verbose, print_all):

    peek_df = pd.read_csv(tabixed_bed_file, sep='\t', nrows=1, header=None, comment='#')
    names = COORDS_COLS5
    if len(peek_df.columns) < len(names):
        msg = f'Invalid blocks file: {tabixed_bed_file}. less than {len(names)} columns.\n'
        msg += f'Run wgbstools convert -L {tabixed_bed_file} -o OUTPUT_REGION_FILE to add the CpG columns'
        raise IllegalArgumentError(msg)

    chroms = GenomeRefPaths().get_chroms()
    params_list = ((tabixed_bed_file, c, pat_file, is_strict, min_len, verbose) for c in chroms)

    p = Pool(num_threads)
    arr = p.starmap(read_blocks_and_test, params_list)
    p.close()
    p.join()

    region_p_val_list = [p_reg for x in arr for p_reg in x]  # flatten
    if not region_p_val_list:
        if verbose:
            eprint(f'[wt bimodal] empty list')
        return
    region_p_val_list = sorted(region_p_val_list, key=lambda elem: elem[1])
    [block_lines, p_vals] = zip(*region_p_val_list)
    accepted_blocks, corrected_p_vals = choose_blocks_by_fdr_bh(p_vals, block_lines, print_all)
    with open(out_file, "w") if out_file != "-" else sys.stdout as f_out:
        for accepted_block, corrected_p_val in zip(accepted_blocks, corrected_p_vals):
            f_out.write(f"{accepted_block}\t{corrected_p_val:,.1e}\n")


def add_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat', help="The input pat file")
    add_GR_args(parser, bed_file=True, required=True)
    add_multi_thread_args(parser)
    parser.add_argument('--strict', action='store_true',
                        help='Truncate reads that start/end outside the given region.')
    parser.add_argument('--min_len', type=int, default=1,
                        help='Only use reads covering at least MIN_LEN CpG sites [1]')
    parser.add_argument('--out_file', '-o', default="-", help="Output file name in which to write results")
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--print_all_regions', action='store_true', help="Print all regions and not only the significant ones.")
    return parser


def parse_args(parser):
    args = parser.parse_args()
    return args


def main():
    """
    Test whether region is bimodal
    """
    parser = add_args()
    args = parse_args(parser)

    if args.bed_file is not None:
        test_multiple_regions(args.bed_file, args.pat, args.threads, args.out_file, args.strict, args.min_len,
                              args.verbose, args.print_all_regions)
    else:
        gr = GenomicRegion(args)
        test_single_region(args.pat, gr.chrom, gr.sites, args.strict, args.min_len)

if __name__ == '__main__':
    main()
