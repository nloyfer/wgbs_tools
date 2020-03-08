#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_single_file, PAT2BETA_TOOL, delete_or_skip, splitextgz, IllegalArgumentError, \
    GenomeRefPaths
import subprocess
import os.path as op
from multiprocessing import Pool
import multiprocessing
from merge import merge_betas

import os


def pat2beta(pat_path, out_dir, args, force=True):
    validate_single_file(pat_path)

    if pat_path.endswith('.pat.gz'):
        cmd = 'gunzip -cd'
    elif pat_path.endswith('.pat'):
        cmd = 'cat'
    else:
        raise IllegalArgumentError('Invalid pat suffix: {}'.format(pat_path))

    out_beta = op.join(out_dir, splitextgz(op.basename(pat_path))[0] + '.beta')
    if not delete_or_skip(out_beta, force):
        return
    nr_sites = GenomeRefPaths(args.genome).nr_sites

    if args.threads > 1 and pat_path.endswith('.pat.gz') and op.isfile(pat_path + '.csi'):
        return mult_pat2beta(pat_path, out_beta, nr_sites, args)

    cmd += ' {} | {} {} {}'.format(pat_path, PAT2BETA_TOOL, out_beta, nr_sites)
    subprocess.check_call(cmd, shell=True)
    return out_beta


def mult_pat2beta(pat_path, out_beta, nr_sites, args):
    processes = []

    with Pool(args.threads) as p:
        chroms = list(GenomeRefPaths(args.genome).get_chrom_cpg_size_table()['chr'])
        for chrom in sorted(chroms):
            beta = '{}.{}.beta'.format(op.splitext(out_beta)[0], chrom)
            params = (chrom, pat_path, beta, nr_sites)
            processes.append(p.apply_async(chr_thread, params))
        p.close()
        p.join()

    # todo: don't save all intermediate beta files. do everything in RAM (edit stdin2beta)
    beta_files = [pr.get() for pr in processes]
    merge_betas(beta_files, out_beta)
    list(map(os.remove, beta_files))
    return out_beta


def chr_thread(chrom, pat, beta, nr_sites):
    cmd = 'tabix {} {} | '.format(pat, chrom)
    cmd += '{} {} {}'.format(PAT2BETA_TOOL, beta, nr_sites)
    subprocess.check_call(cmd, shell=True, stderr=subprocess.PIPE)
    return beta


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat_path', help='A pat[.gz] file')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if existed')
    parser.add_argument('-o', '--out_dir', help='Output directory for the beta file. [.]', default='.')
    parser.add_argument('--genome', help='Genome reference name. Default is hg19.', default='hg19')
    parser.add_argument('-@', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='Number of threads to use (default: multiprocessing.cpu_count)')
    return parser.parse_args()


def main():
    """
    Generate a beta file from a pat file
    """
    args = parse_args()
    pat2beta(args.pat_path, args.out_dir, args, args.force)


if __name__ == '__main__':
    main()
