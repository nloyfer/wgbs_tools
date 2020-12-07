import numpy as np
import pandas as pd
import argparse
from utils_wgbs import add_GR_args
from genomic_region import GenomicRegion
from segmentor_client import insert_genomic_loci, dump_table


def process_block_file(in_file, args, out_file):
    cpg_starts = []
    cpg_ends = []
    with open(in_file, 'rt') as block_lines:
        for block_line in block_lines:
            block_tokens = block_line.split("\t")
            cpg_start = block_tokens[1]
            cpg_end = block_tokens[2]
            cpg_starts.append(cpg_start)
            cpg_ends.append(cpg_end)
    start_arr = np.array([int(start_ind) for start_ind in cpg_starts])
    end_arr = np.array([int(end_ind) for end_ind in cpg_ends])
    df = pd.DataFrame({'startCpG': start_arr, 'endCpG': end_arr})
    gr = GenomicRegion(args)
    df = insert_genomic_loci(df, gr)
    dump_table(df, out_file)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file')
    parser.add_argument('--out_file')
    add_GR_args(parser)
    return parser.parse_args()


def main():
    args = parse_args()
    process_block_file(args.in_file, args, args.out_file)
    from index_wgbs import Indxer
    Indxer(args.out_file).run()


if __name__ == '__main__':
    main()