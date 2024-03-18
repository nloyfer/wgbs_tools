#!/usr/bin/python3 -u

import argparse
import subprocess
import os.path as op
from multiprocessing import Pool
import numpy as np
from utils_wgbs import validate_single_file, pat2beta_tool, \
        delete_or_skip, splitextgz, IllegalArgumentError, \
        GenomeRefPaths, trim_to_uint8, add_multi_thread_args, \
        validate_local_exe


def pat2beta(pat_path, out_dir, args, force=True):
    validate_single_file(pat_path)

    if pat_path.endswith('.pat.gz'):
        cmd = 'gunzip -cd'
    elif pat_path.endswith('.pat'):
        cmd = 'cat'
    else:
        raise IllegalArgumentError(f'Invalid pat suffix: {pat_path}')

    suff = '.lbeta' if args.lbeta else '.beta'
    out_beta = op.join(out_dir, splitextgz(op.basename(pat_path))[0] + suff)
    if not delete_or_skip(out_beta, force):
        return

    if args.threads > 1 and pat_path.endswith('.pat.gz') and op.isfile(pat_path + '.csi'):
        arr = mult_pat2beta(pat_path, args)
    else:
        nr_sites = GenomeRefPaths(args.genome).get_nr_sites()
        cmd += f' {pat_path} | {pat2beta_tool} {1} {nr_sites + 1}'
        x = subprocess.check_output(cmd, shell=True).decode()
        arr = np.fromstring(x, dtype=int, sep=' ').reshape((-1, 2))

    trim_to_uint8(arr, args.lbeta).tofile(out_beta)
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
    cmd += f'{pat2beta_tool} {start} {end}'
    x = subprocess.check_output(cmd, shell=True).decode()
    x = np.fromstring(x, dtype=int, sep=' ').reshape((-1, 2))
    return x


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat_paths', help='pat[.gz] files', nargs='+')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if existed')
    parser.add_argument('-o', '--out_dir', help='Output directory for the beta file. [.]', default='.')
    parser.add_argument('-l', '--lbeta', action='store_true', help='Use lbeta file (uint16) instead of beta (uint8)')
    parser.add_argument('--genome', help='Genome reference name.')
    add_multi_thread_args(parser)
    return parser.parse_args()


def main():
    """
    Generate a beta file from a pat file
    """
    args = parse_args()
    validate_local_exe(pat2beta_tool)
    for pat in args.pat_paths:
        pat2beta(pat, args.out_dir, args, args.force)


if __name__ == '__main__':
    main()
