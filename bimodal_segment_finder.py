import argparse
import sys
import numpy as np
import pandas as pd
import math
from multiprocessing import Pool
import multiprocessing
import subprocess


def get_beta_array(file_name):
    return np.fromfile(file_name, dtype=np.uint8).reshape((-1, 2))


def get_dis_from_half_list(blocks_df, beta_arrays):
    avg_dis_from_half_list = []
    for index, row in blocks_df.iterrows():
        start_ind = row.cpg_start - 1
        end_ind = row.cpg_end - 1
        avg_dist_from_half = get_avg_dist_from_half(start_ind, end_ind, beta_arrays)
        avg_dis_from_half_list.append(avg_dist_from_half)
    return avg_dis_from_half_list


def get_avg_dist_from_half(start_ind, end_ind, beta_arrays, coverage):
    average_dist_from_half = 0
    total = 0
    cur_coverage = 0
    for b_array in beta_arrays:
        cur_b_array = b_array[start_ind:end_ind, :]
        for row in cur_b_array:
            methyl_prop = 0 if row[1] == 0 else row[0] / row[1]
            average_dist_from_half += abs(methyl_prop - 0.5)
            cur_coverage += row[1]
            total += 1

    return 1 if cur_coverage <= coverage else average_dist_from_half / total


def scan_blocks_for_switches(block_file, beta_file, out_path, distance_threshold):
    blocks_df = pd.read_csv(block_file, sep='\t', names=['chrom', 'start', 'end', 'cpg_start', 'cpg_end'])
    b_array = [get_beta_array(beta_file)]
    dist_from_half_list = get_dis_from_half_list(blocks_df, b_array)
    blocks_df['dist_from_half'] = dist_from_half_list
    new_blocks = blocks_df[blocks_df['dist_from_half'] < distance_threshold]
    new_blocks.to_csv(out_path, sep='\t', index=False, header=None)


def scan_block_for_bimodal_start(row, beta_arrays, distance_threshold, window_size=10, consecutive_wins=10, coverage=8):
    start_ind = row.cpg_start - 1
    end_ind = row.cpg_end + 1
    cur_ind = start_ind
    num_wins = 0
    real_start = None
    while cur_ind + window_size <= end_ind:
        average_dist_from_half = get_avg_dist_from_half(cur_ind, cur_ind + window_size, beta_arrays, coverage)
        if average_dist_from_half < distance_threshold:
            num_wins += 1
        else:
            num_wins = 0
        if num_wins >= consecutive_wins:
            real_start = cur_ind - consecutive_wins
            break
        cur_ind += 1
    if real_start:
        real_start += 1
    else:
        real_start = math.nan
    return real_start


def scan_block_for_bimodal_end(row, beta_arrays, distance_threshold, window_size=10, consecutive_wins=10, coverage=8):
    start_ind = row.cpg_start - 1
    end_ind = row.cpg_end + 1
    cur_ind = end_ind
    num_wins = 0
    real_end = None
    while cur_ind - window_size >= start_ind:
        average_dist_from_half = get_avg_dist_from_half(cur_ind - window_size, cur_ind, beta_arrays, coverage)
        if average_dist_from_half < distance_threshold:
            num_wins += 1
        else:
            num_wins = 0
        if num_wins >= consecutive_wins:
            real_end = cur_ind + consecutive_wins
            break
        cur_ind -= 1
    if real_end:
        real_end += 1
    else:
        real_end = math.nan
    return real_end


def check_block_min_u_m_proportion(row, pat_file_name, num_sites_per_read=4, homog_cutoff=0.75, min_u_m_proportion=0.3):
    start_ind = row.real_start
    end_ind = row.real_end
    chrom = row.chrom
    tabix_cmd = "tabix {} {}:{}-{}".format(pat_file_name, chrom, start_ind, end_ind)
    p = subprocess.Popen(tabix_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    u_count = 0
    m_count = 0
    x_count = 0
    for line in iter(p.stdout.readline, ''):
        line = line.decode()
        line = line.rstrip()
        if len(line) == 0:
            break
        tokens = line.split("\t")
        cur_count = int(tokens[3])
        pattern = tokens[2]
        cur_u_count = 0
        cur_m_count = 0
        for el in pattern:
            if el == 'C':
                cur_m_count += 1
            elif el == 'T':
                cur_u_count += 1
        total = cur_m_count + cur_u_count
        u_proportion = 0 if total == 0 else cur_u_count / total
        m_proportion = 0 if total == 0 else cur_m_count / total
        if total >= num_sites_per_read:
            if u_proportion >= homog_cutoff:
                u_count += cur_count
            elif m_proportion >= homog_cutoff:
                m_count += cur_count
            else:
                x_count += cur_count
    total_reads = u_count + m_count
    reads_u_proportion = 0 if total_reads == 0 else u_count / total_reads
    reads_m_proportion = 0 if total_reads == 0 else m_count / total_reads
    if min(reads_m_proportion, reads_u_proportion) >= min_u_m_proportion:
        return 1
    else:
        return 0


def run_apply_start_end(df, col_name, method_to_run, beta_arrays, distance_threshold, window_size, consecutive_wins,
                        coverage):
    x = np.array(df.apply(
        lambda row: method_to_run(row, beta_arrays, distance_threshold, window_size, consecutive_wins,
                                  coverage), axis=1))
    df[col_name] = x
    df = df.dropna()
    return df


def run_apply_check_block(df, pat_file):
    df['above_u_m_threshold'] = df.apply(
        lambda row: check_block_min_u_m_proportion(row, pat_file), axis=1)
    df = df[df.above_u_m_threshold == 1]
    df = df.drop('above_u_m_threshold', 1)
    return df


def run_apply(df, col_name):
    df[col_name] = df.apply(lambda row: 1, axis=1)
    df = df.dropna()
    return df


def process_df(df_file):
    blocks_df = pd.read_csv(df_file, sep='\t', names=['a', 'b', 'c', 'd', 'e'])
    group_size = 2000
    df_list = [blocks_df.iloc[i:i + group_size] for i in range(0, len(blocks_df), group_size)]
    params_1 = [(df, 'new_a') for df in df_list]
    p = Pool(1)
    arr_from_start = p.starmap(run_apply, params_1)
    p.close()
    p.join()
    params_2 = [(df, 'new_b') for df in arr_from_start]
    p = Pool(1)
    arr_from_start = p.starmap(run_apply, params_2)
    p.close()
    p.join()
    blocks_df = pd.concat(arr_from_start)
    print(blocks_df.shape)


def iterate_on_blocks_for_bimodal_region(block_file, beta_file, out_path, pat_file,
                                         distance_threshold, num_threads, window_size=10,
                                         consecutive_wins=10, coverage=8):
    blocks_df = pd.read_csv(block_file, sep='\t', names=['chrom', 'start', 'end', 'cpg_start', 'cpg_end'])
    beta_arrays = [get_beta_array(beta_file)]
    group_size = 2000
    df_list = [blocks_df.iloc[i:i + group_size] for i in range(0, len(blocks_df), group_size)]
    params_for_scan_from_start = [(df, 'real_start', scan_block_for_bimodal_start, beta_arrays, distance_threshold,
                                   window_size, consecutive_wins, coverage) for df in df_list]
    print("starting to scan for start of bimodal window")
    p = Pool(num_threads)
    arr_from_start = p.starmap(run_apply_start_end, params_for_scan_from_start)
    p.close()
    p.join()

    params_for_scan_from_end = [
        (df, 'real_end', scan_block_for_bimodal_end, beta_arrays, distance_threshold, window_size, consecutive_wins,
         coverage) for
        df in arr_from_start]
    print("starting to scan for end of bimodal window")
    p = Pool(num_threads)
    arr_from_end = p.starmap(run_apply_start_end, params_for_scan_from_end)
    p.close()
    p.join()
    # blocks_df['real_start'] = blocks_df.apply(
    #     lambda row: scan_block_for_bimodal_start(row, beta_arrays, distance_threshold, window_size, consecutive_wins,
    #                                              coverage), axis=1)
    # blocks_df['real_end'] = blocks_df.apply(
    #     lambda row: scan_block_for_bimodal_end(row, beta_arrays, distance_threshold, window_size, consecutive_wins,
    #                                            coverage), axis=1)
    # blocks_df = blocks_df.dropna()
    params_for_check_block = [
        (df, pat_file) for
        df in arr_from_end]
    # blocks_df['above_u_m_threshold'] = blocks_df.apply(lambda row: check_block_min_u_m_proportion(row, pat_file),
    #                                                    axis=1)
    # blocks_df = blocks_df[blocks_df.above_u_m_threshold == 1]
    # blocks_df = blocks_df.drop('above_u_m_threshold', 1)
    print("ruling out non-bimodal windows")
    p = Pool(num_threads)
    arr_final_blocks = p.starmap(run_apply_check_block, params_for_check_block)
    p.close()
    p.join()
    blocks_df = pd.concat(arr_final_blocks)
    blocks_df.astype({'real_start': 'int32', 'real_end': 'int32'}).to_csv(out_path, sep='\t', index=False, header=None)


def parse_args():
    parser = argparse.ArgumentParser()
    # betas_or_file = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('--beta_file', default="/cs/cbio/jon/grail_atlas/data/Prostate-Epithelial-Z000000RV.beta",
                        help="Input beta file")
    parser.add_argument('--block_file', default="/cs/cbio/jon/segmentation_files/Prostate-Epithelial-Z000000RV.blocks.tsv",
                        help="Input block file")
    parser.add_argument('--pat_file', default="/cs/cbio/jon/grail_atlas/data/Prostate-Epithelial-Z000000RV.pat.gz",
                        help="Input pat file")
    parser.add_argument('--dist_threshold', default=0.22, help="Threshold on the average distance from 0.5 beta value "
                                                               "per block. Default is 0.22")
    parser.add_argument('--min_num_cpgs', default=8, help="Minimum number of observed CpG sites to include block. Default is 8.")
    parser.add_argument('-w', '--window_size', default=6, help="window size. Default is 6.")
    parser.add_argument('-@', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='Number of threads to use (default: multiprocessing.cpu_count)')
    parser.add_argument('-o', '--out_path', default="/cs/cbio/jon/segmentation_files/Prostate-Epithelial-Z000000RV.blocks.bimodal.tsv",
                        # sys.stdout,  # "/cs/cbio/jon/trial_blocks.tsv",
                        help='output path [stdout]')
    return parser.parse_args()


def main():
    args = parse_args()
    iterate_on_blocks_for_bimodal_region(args.block_file, args.beta_file, args.out_path, args.pat_file,
                                         float(args.dist_threshold), args.threads, window_size=int(args.window_size),
                                         consecutive_wins=int(args.window_size),
                                         coverage=int(args.min_num_cpgs))
    # process_df(args.block_file)


if __name__ == '__main__':
    main()
