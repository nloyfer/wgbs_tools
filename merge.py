#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_files_list, splitextgz, delete_or_skip, trim_to_uint8, load_beta_data, \
    collapse_pat_script, IllegalArgumentError
import subprocess
import numpy as np
import os.path as op
import sys
import os
from index_wgbs import index_single_file


def merge_sort_open_pats(open_pats, merged_path, remove_open_pats=True):
    """
    Merge-sort a list of pat files into one pat. Then bgzip and index it.
    Remove the open_pats
    :param open_pats: list of sorted pat files, not compressed.
    :param merged_path: final output pat path, not compressed (ends with 'pat')
    :return: path of compressed merged pat
    """
    # merge-sort all pats:
    pat_with_dups = merged_path + '.with_dups'          # i.e "prefix.pat.with_dups"
    cmd = 'sort -m -k2,2n -k3,3 ' + ' '.join(open_pats)
    cmd += ' -o {} '.format(pat_with_dups)
    if subprocess.call(cmd, shell=True):
        raise IllegalArgumentError('Failed in merge sort')

    # remove the open_pats files
    if remove_open_pats:
        for tmp in open_pats:
            os.remove(tmp)

    # collapse "almost duplicated" lines:
    if subprocess.call('{} {} > {}'.format(collapse_pat_script, pat_with_dups, merged_path), shell=True):
        raise IllegalArgumentError('Failed collapsing pat')

    os.remove(pat_with_dups)

    # bgzip and index the output
    if index_single_file(merged_path, force=True):
        raise IllegalArgumentError('Failed indexing pat')

    return merged_path + '.gz'


def merge_pats(pats, outpat):
    """
    Merge all pats as follows:
        - decompress all pats to temporary files.
        - Merge-sort them to a united uncompressed pat file with "almost duplicated" lines
        - collapse "almost duplicated" lines
        - bgzip and index result
    :param pats: A list of sorted .pat.gz files
    :param outpat: the output path. ends with ".pat.gz"
    :return: 0 iff success
    """
    # gunzip all parts:
    tmps = []
    for pat in pats:
        pat_tmp = op.splitext(pat)[0] + '.tmp'      # i.e "part.pat.tmp"
        tmps.append(pat_tmp)
        subprocess.call('gunzip -cd {} > {}'.format(pat, pat_tmp), shell=True)

    # merge-sort all pats:
    outpat_nogzip = op.splitext(outpat)[0]          # i.e "prefix.pat"
    merge_sort_open_pats(tmps, outpat_nogzip)


def merge_betas(betas, opath):
    """
    Merge all betas by summing their values element-wise, while keeping the dimensions
    :param betas: list of beta files
    :param opath: merged beta file
    """
    data = load_beta_data(betas[0]).astype(np.int)
    for b in betas[1:]:
        data += load_beta_data(b)

    # Trim / normalize to range [0, 256)
    data = trim_to_uint8(data)
    # Dump
    data.tofile(opath)
    return data


def merge_unqs(unqs, opath):
    raise NotImplemented


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+')
    parser.add_argument('-p', '--prefix', help='Prefix of output file', required=True)
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if existed')
    # parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    return args


def main():
    """
    Merge files.
    Accumulate all reads / observations from multiple (>=2) input files,
    and output a single file of the same format.
    Supported formats: pat.gz, beta
    """
    args = parse_args()

    # validate input files
    input_files = args.input_files
    validate_files_list(input_files, min_len=2)

    # construct output path
    out_path = args.prefix + splitextgz(args.input_files[0])[1]
    if not delete_or_skip(out_path, args.force):
        return

    params = (input_files, out_path)

    if input_files[0].endswith('.beta'):
        merge_betas(*params)
    elif input_files[0].endswith('.pat.gz'):
        merge_pats(*params)
    elif input_files[0].endswith('.unq.gz'):
        merge_unqs(*params)
    else:
        print('Unknown input format:', input_files[0])
        return


if __name__ == '__main__':
    main()
