#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_files_list, splitextgz, delete_or_skip, trim_to_uint8, load_beta_data, \
    collapse_pat_script, IllegalArgumentError
import subprocess
import numpy as np
import os.path as op
import sys
import os
from index_wgbs import Indxer
import pandas as pd


class MergePats:
    def __init__(self, pats, out_nogzip, labels):
        self.pats = pats
        self.out_nogzip = out_nogzip
        self.labels = labels

    def merge_pats(self, fast=False):
        self.fast_merge_pats(None) if fast else self.slow_merge_pats()

    def slow_merge_pats(self):
        """
        Merge all pats as follows:
            - decompress all pats to temporary files.
            - Merge-sort them to a united uncompressed pat file with "almost duplicated" lines
            - collapse "almost duplicated" lines
            - bgzip and index result
        :return: 0 iff success
        """

        # gunzip all parts:
        tmps = []
        for pat in self.pats:
            pat_tmp = op.splitext(pat)[0] + '.tmp'  # i.e "part.pat.tmp"
            tmps.append(pat_tmp)
            subprocess.call('gunzip -cd {} > {}'.format(pat, pat_tmp), shell=True)

        # merge-sort all pats:
        merge_sort_open_pats(tmps, self.out_nogzip)

    def compose_view_cmd(self, i, view_flags, labels):
        pat = self.pats[i]
        view_cmd = ''
        if view_flags is None:
            view_cmd += ' <(zcat {p}'.format(p=pat)
        else:
            view_cmd += ' <(wgbs_tools view {flags} {p} '.format(p=pat,
                                                                 flags=view_flags[i])
        tagcmd = ''
        if labels is not None:
            tagcmd = ' | sed s/\$/\'\t\'{}/'.format(labels[i])
        view_cmd += '{})'.format(tagcmd)
        return view_cmd

    def fast_merge_pats(self, view_flags=None):
        """
        Use bash script, not intermediate files
        """
        cmd = 'sort -m -k2,2n -k3,3'
        for i in range(len(self.pats)):
            cmd += self.compose_view_cmd(i, view_flags, self.labels)
        cmd += ' | {} - '.format(collapse_pat_script)
        cmd += ' | bgzip > {}.gz'.format(self.out_nogzip)
        cmd = '/bin/bash -c "{}"'.format(cmd)
        subprocess.check_call(cmd, shell=True)
        if not op.isfile(self.out_nogzip + '.gz'):
            raise IllegalArgumentError('Error: failed to create file {}.gz'.format(self.out_nogzip))
        Indxer(self.out_nogzip + '.gz', force=True).run()


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
    subprocess.check_call(cmd, shell=True)

    # remove the open_pats files
    if remove_open_pats:
        for tmp in open_pats:
            if op.isfile(tmp):
                os.remove(tmp)

    # collapse "almost duplicated" lines:   # todo: use python instead of perl. pandas duplicates, then loop.
    subprocess.check_call('{} {} > {}'.format(collapse_pat_script, pat_with_dups, merged_path), shell=True)

    #df = pd.read_csv(pat_with_dups, sep='\t', header=None)
    #collapse_adup_df(df)
    os.remove(pat_with_dups)

    # bgzip and index the output
    Indxer(merged_path, force=True).run()

    return merged_path + '.gz'






#def compose_cols_names(pat_path, names):
#    peek = pd.read_csv(pat_path, sep='\t', header=None, nrows=1)
#    i = 0
#    while len(peek.columns) > len(names):
#        names += ['tag{}'.format(i)]
#        i += 1
#    return names


def collapse_adup_df(df):
    """ given a dataframe, collapse almost duplicated lines and return it """#    
    print(df)
    names_to_test = list(df.columns).remove(2)
    bs = df.duplicates(subset=names_to_test)
    print(bs)
    print(df)
    return df

#def collapse_adup(pat_in, pat_out):
#    for cf in df_generator(pat_in, 10):
#        collapse_adup_df(cf)
#        cf.to_csv(pat_out, sep='\t', header=None, index=None, mode='a')
#    # todo: bgzip and index
#    return


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


# class MergeUnqs:
#     def __init__(self, unqs, out_nogzip):
#         self.unqs = unqs
#         self.out_nogzip = out_nogzip
#
#     def compose_view_cmd(self, i, view_flags):
#         if view_flags is None:
#             return ' <(gunzip -cd {p})'.format(p=self.unqs[i])
#         else:
#             return ' <(wgbs_tools view {flags} {p})'.format(p=self.unqs[i], flags=view_flags[i])
#
#     def fast_merge_unqs(self, view_flags=None):
#         """
#         Use bash script, with no intermediate files
#         """
#         cmd = 'sort -m -k2,2n -k3,3'
#         for i in range(len(self.unqs)):
#             cmd += self.compose_view_cmd(i, view_flags)
#         cmd +=  ' | {} - '.format(collapse_pat_script)
#         # cmd += ' | bgzip > {o}.gz && tabix -C -b 2 -e 2 {o}.gz'.format(o=self.out_nogzip)
#         cmd += ' | bgzip > {}.gz'.format(self.out_nogzip)
#         cmd = '/bin/bash -c "{}"'.format(cmd)
#         # print(cmd)
#         subprocess.check_call(cmd, shell=True)
#         Indxer(self.out_nogzip + '.gz', force=True).run()

def merge_unqs():
    raise NotImplemented


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+')
    parser.add_argument('-p', '--prefix', help='Prefix of output file', required=True)
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if existed')
    parser.add_argument('--fast', action='store_true', help='pat: Use bash for better performance')
    # parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--labels', nargs='+', help='labels for the mixed reads. '
                                                    'Default is None')
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

    files_type = splitextgz(input_files[0])[1][1:]
    labels = None if not args.labels else args.labels

    if files_type in ('beta', 'bin'):
        merge_betas(input_files, out_path)
    elif files_type == 'pat.gz':
        MergePats(input_files, args.prefix + '.pat', labels).merge_pats(args.fast)
    elif files_type == 'unq.gz':
        merge_unqs()
    else:
        print('Unknown input format:', input_files[0])
        return


if __name__ == '__main__':
    main()
