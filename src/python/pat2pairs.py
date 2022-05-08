#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_single_file,  delete_or_skip, splitextgz, IllegalArgumentError, \
    GenomeRefPaths, add_multi_thread_args
import subprocess
import os.path as op
from multiprocessing import Pool
import multiprocessing
import numpy as np
import os
from pathlib import Path

path = Path(op.realpath(__file__))
DIR = str(path.parent)
SRC_DIR = op.join(path.parent.parent.parent, 'src/')
PAT2PAIR_TOOL = SRC_DIR + 'pat2beta/stdin2pairs'

def trim_to_uint8(data):
    # Trim / normalize to range [0, 256)

    max_val = 2 ** 8 - 1
    big_inds = np.argwhere(data.max(axis=1) > max_val).flatten()
    data[big_inds, :] = data[big_inds, :] / data.max(axis=1)[big_inds][:, None] * max_val
    return data.astype(np.uint8)

def pat2pairs(pat_path, out_dir, args, force=True):
    validate_single_file(pat_path)

    if pat_path.endswith('.pat.gz'):
        cmd = 'gunzip -cd'
    elif pat_path.endswith('.pat'):
        cmd = 'cat'
    else:
        raise IllegalArgumentError(f'Invalid pat suffix: {pat_path}')

    suff = '.pairs'
    out_beta = op.join(out_dir, splitextgz(op.basename(pat_path))[0] + suff)
    if not delete_or_skip(out_beta, force):
        return

    if args.threads > 1 and pat_path.endswith('.pat.gz') and op.isfile(pat_path + '.csi'):
        arr = mult_pat2beta(pat_path, args)
    else:
        nr_sites = GenomeRefPaths(args.genome).nr_sites
        cmd += f' {pat_path} | {PAT2PAIR_TOOL} {1} {nr_sites + 1}'
        x = subprocess.check_output(cmd, shell=True).decode()
        arr = np.fromstring(x, dtype=int, sep=' ').reshape((-1, 4))

    trim_to_uint8(arr).tofile(out_beta)
    return out_beta


def mult_pat2beta(pat_path, args):
    processes = []
    with Pool(args.threads) as p:
        ct = GenomeRefPaths(args.genome).get_chrom_cpg_size_table()
        x = np.cumsum([0] + list(ct['size'])) + 1
        chroms = list(ct['chr'])
        for i, chrom in enumerate(chroms):
            start = x[i]
            end = x[i + 1]
            params = (pat_path, chrom, start, end)
            processes.append(p.apply_async(chr_thread, params))
        p.close()
        p.join()

    beta_files = [pr.get() for pr in processes]
    res = np.concatenate(beta_files, axis = 0)
    return res


def chr_thread(pat, chrom, start, end):
    cmd = f'tabix {pat} {chrom} | '
    cmd += f'{PAT2PAIR_TOOL} {start} {end}'
    x = subprocess.check_output(cmd, shell=True).decode()
    x = np.fromstring(x, dtype=int, sep=' ').reshape((-1, 4))
    return x


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat_paths', help='pat[.gz] files', nargs='+')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if existed')
    parser.add_argument('-o', '--out_dir', help='Output directory for the beta file. [.]', default='.')
    parser.add_argument('--genome', help='Genome reference name. Default is hg19.', default='hg19')
    add_multi_thread_args(parser)
    return parser.parse_args()


def main():
    """
    Generate a beta file from a pat file
    """
    args = parse_args()
    for pat in args.pat_paths:
        pat2pairs(pat, args.out_dir, args, args.force)


if __name__ == '__main__':
    main()
