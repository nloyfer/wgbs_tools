#!/usr/bin/python3 -u

import argparse
import os
import numpy as np
import os.path as op
import pandas as pd
from multiprocessing import Pool
import sys
import subprocess
from io import StringIO
from beta_to_blocks import load_blocks_file, is_block_file_nice
from utils_wgbs import IllegalArgumentError, add_multi_thread_args, \
        homog_tool, main_script, splitextgz, GenomeRefPaths, validate_file_list \
        COORDS_COLS5


def homog_log(*args, **kwargs):
    print('[ wt homog ]', *args, file=sys.stderr, **kwargs)


import pdb
# ForkedPdb().set_trace()
class ForkedPdb(pdb.Pdb):
    """A Pdb subclass that may be used
    from a forked multiprocessing child

    """
    def interaction(self, *args, **kwargs):
        _stdin = sys.stdin
        try:
            sys.stdin = open('/dev/stdin')
            pdb.Pdb.interaction(self, *args, **kwargs)
        finally:
            sys.stdin = _stdin


######################################################
#                                                    #
#    wrap the c++ tool                               #
#                                                    #
######################################################


def trim_uxm_to_uint8(data, nr_bits):
    data = data.copy()
    if nr_bits == 16:
        dtype = np.uint16
    else:
        dtype = np.uint8
    max_val = 2 ** nr_bits - 1
    big_inds = np.argwhere(data.max(axis=1) > max_val).flatten()
    data[big_inds, :] = data[big_inds, :] / data.max(axis=1)[big_inds][:, None] * max_val
    res = data.astype(dtype)
    return res

def homog_chrom(pat, name, blocks, blocks_path, chrom, rates_cmd, view_full):
    if view_full:
        cmd = f'tabix {pat} {chrom}'
    else:
        cmd = f'{main_script} cview {pat} -L {blocks_path}'
    cmd += f' | {homog_tool} -b {blocks_path} -n {name}.{chrom} {rates_cmd}'
    cmd += f' --chrom {chrom}'
    txt = subprocess.check_output(cmd, shell=True).decode()
    names = ['U', 'X', 'M']
    df = pd.read_csv(StringIO(txt), sep='\t', header=None, names=names)
    if df.values.sum() == 0:
        homog_log(f' [ {name}.{chrom} ] WARNING: all zeros!')

    df = pd.concat([blocks.reset_index(drop=True), df], axis=1)
    return df

def should_be_skipped(force, bin_path, bed_path, binary, bed):
    if force:
        return False
    if (not op.isfile(bin_path)) and (not op.isfile(bed_path)):
        return False
    rerun_bed = (not op.isfile(bed_path)) and bed
    rerun_bin = (not op.isfile(bin_path)) and binary
    return not (rerun_bin or rerun_bed)

def homog_process(pat, blocks, args):
    name = splitextgz(op.basename(pat))[0]
    prefix = op.join(args.out_dir, name)
    bin_path = prefix + '.uxm'
    bed_path = prefix + '.uxm.bed.gz'
    bed = args.bed
    binary = args.binary or (not args.binary and not bed)
    if should_be_skipped(args.force, bin_path, bed_path, binary, bed):
        homog_log(f'skipping {name}. Use -f to overwrite')
        return

    homog_log(f' [ {name} ] starting')

    # generate rate_cmd:
    l = args.rlen
    rate_cmd = f' -l {l} -r '
    if args.thresholds:
        rate_cmd += f'0,{args.thresholds},1'
    else:
        th1 = round(1 - (l - 1) / l, 3) + 0.001
        th2 = round((l - 1) / l, 3)
        rate_cmd += f'0,{th1},{th2},1 '

    # for a long marker file (>10K marker), 
    # parse the whole pat file instead of running "cview -L BED"
    view_full = blocks.shape[0] > 1e4

    cf = GenomeRefPaths(args.genome).get_chrom_cpg_size_table()
    chroms = sorted(set(cf.chr) & set(blocks.chr))
    p = Pool(args.threads)
    params = [(pat, name, blocks[blocks['chr'] ==c],
               args.blocks_file, c, rate_cmd,
               view_full) for c in chroms]
    arr = p.starmap(homog_chrom, params)
    p.close()
    p.join()

    df = pd.concat(arr).reset_index(drop=True)
    df = blocks.merge(df, how='left', on=COORDS_COLS5)

    if binary:
        trim_uxm_to_uint8(df[list('UXM')].values, args.nr_bits).tofile(bin_path)
    if bed:
        df.to_csv(bed_path, sep='\t', header=None, index=None)
    return df


def main():  # TODO: this is 8x slower than simply running the CPP tool on a single thread. Fix this.
    """
    Generage homog files. Given a blocks file and pat[s],
    count the number of U,X,M reads for each block for each file
    """

    args = parse_args()
    if args.nr_bits not in (8 , 16):
        raise IllegalArgumentError('nr_bits must be in {8, 16}')
    if args.rlen < 3:
        raise IllegalArgumentError('rlen must be >= 3')
    if args.thresholds is not None:
        th = args.thresholds.split(',')
        if not len(th) == 2: # and th[0].is_number():
            raise IllegalArgumentError('Invalid thresholds')
        th = float(th[0]), float(th[1])
        if not (1 > th[1] > th[0] > 0):
            raise IllegalArgumentError('Invalid thresholds')
    pats = args.input_files
    validate_file_list(pats, '.pat.gz')

    # load blocks:
    blocks_df = load_blocks_file(args.blocks_file)
    is_nice, msg = is_block_file_nice(blocks_df)
    if not is_nice:
        homog_log(msg)
        raise IllegalArgumentError(f'Invalid blocks file: {args.blocks_file}')

    for pat in sorted(pats):
        homog_process(pat, blocks_df, args)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+', help='one or more pat files')
    parser.add_argument('-b', '--blocks_file', help='blocks path', required=True)
    parser.add_argument('-o', '--out_dir', help='output directory. Default is "."', default='.')
    parser.add_argument('--force', '-f', action='store_true', help='Overwrite files if exist')
    parser.add_argument('--binary', action='store_true', help='Output binary files (uint8)')
    parser.add_argument('--bed', action='store_true', help='Output bed file')
    parser.add_argument('--genome', help='Genome reference name.')
    parser.add_argument('--nr_bits', type=int, default=8,
            help='For binary output, specify number of bits for the output format - 8 or 16. ' \
                 '(e.g. 8 statnds for uint8, which means values are trimmed to [0, 255])')
    parser.add_argument('--thresholds', '-t',
            help='UXM thresholds, LOW,HIGH. E.g, "0.3334,0.666".\n')
    parser.add_argument('--rlen', '-l', type=int, default=3,
            help='Minimal read length (in CpGs) to consider. Default is 3')
    parser.add_argument('--debug', '-d', action='store_true')
    add_multi_thread_args(parser)

    return parser.parse_args()


if __name__ == '__main__':
    main()
