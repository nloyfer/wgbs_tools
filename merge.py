#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_files_list, splitextgz, delete_or_skip, trim_to_uint8, load_beta_data, \
    collapse_pat_script, IllegalArgumentError, main_script, eprint
import subprocess
import numpy as np
import os.path as op
import sys
import os
from index_wgbs import Indxer


class MergePats:
    def __init__(self, pats, out_nogzip, labels, args):
        self.args = args
        self.pats = pats
        self.out_nogzip = out_nogzip
        self.labels = labels

    def merge_pats(self):
        self.fast_merge_pats()

    def compose_view_cmd(self, i, view_flags):
        if view_flags is None:
            view_cmd = ' <(gunzip -c'
        else:
            view_cmd = ' <({wt} view {flags}'.format(wt=main_script, flags=view_flags[i])
        view_cmd += ' {}'.format(self.pats[i])
        tagcmd = ''
        if self.labels is not None:
            tagcmd = ' | sed s/\$/\'\t\'{}/'.format(self.labels[i])
        view_cmd += '{})'.format(tagcmd)
        return view_cmd

    def fast_merge_pats(self, view_flags=None):
        """ Use piping and sort -m to merge pat files w/o intermediate files """
        cmd = 'sort -m -k2,2n -k3,3'
        if self.args.temp_dir:
            temp_dir = self.args.temp_dir
            if not op.isdir(temp_dir):
                eprint('Invalid temp dir: {}. Ignoring it'.format(temp_dir))
            else:
                cmd += ' -T {} '.format(temp_dir)
        for i in range(len(self.pats)):
            cmd += self.compose_view_cmd(i, view_flags)
        cmd += ' | {} - '.format(collapse_pat_script)
        cmd += ' | bgzip > {}.gz'.format(self.out_nogzip)
        cmd = '/bin/bash -c "{}"'.format(cmd)
        subprocess.check_call(cmd, shell=True)

        if not op.isfile(self.out_nogzip + '.gz'):
            raise IllegalArgumentError('Error: failed to create file {}.gz'.format(self.out_nogzip))
        Indxer(self.out_nogzip + '.gz', force=self.args.force).run()


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


def merge_unqs():
    raise NotImplemented


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+')
    parser.add_argument('-p', '--prefix', help='Prefix of output file', required=True)
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if existed')
    parser.add_argument('-T', '--temp_dir', help='passed to "sort -m". Useful for merging very large pat files')
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

    if op.realpath(out_path) in [op.realpath(p) for p in args.input_files]:
        eprint('[merge] Error output path is identical ' \
                'to one of the input files {}'.format(out_path))
        return

    if not delete_or_skip(out_path, args.force):
        return

    files_type = splitextgz(input_files[0])[1][1:]

    if files_type in ('beta', 'bin'):
        merge_betas(input_files, out_path)
    elif files_type == 'pat.gz':
        MergePats(input_files, args.prefix + '.pat', args.labels, args).merge_pats()
    elif files_type == 'unq.gz':
        merge_unqs()
    else:
        print('Unknown input format:', input_files[0])
        return


if __name__ == '__main__':
    main()
