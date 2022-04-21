#!/usr/bin/python3 -u

import argparse
import subprocess
import numpy as np
import os.path as op
import sys
import os
from index import Indxer
from utils_wgbs import validate_file_list, splitextgz, delete_or_skip, trim_to_uint8, load_beta_data, \
        collapse_pat_script, IllegalArgumentError, main_script, eprint, add_GR_args
from genomic_region import GenomicRegion


class MergePats:
    def __init__(self, pats, outpath, labels, args):
        self.args = args
        self.pats = pats
        validate_file_list(self.pats, force_suff='.pat.gz')
        self.outpath = outpath
        self.labels = labels

    def merge_pats(self):
        view_flags = []
        for i in range(len(self.pats)):
            v = ' '
            if self.args.strict:
                v += ' --strict'
            if self.args.min_len:
                v += ' --min_len {}'.format(self.args.min_len)
            if self.args.bed_file is not None:
                v += ' -L {}'.format(self.args.bed_file)
            gr = GenomicRegion(self.args)
            if not gr.is_whole():
                v += ' -s {}-{}'.format(*gr.sites)
            # v += ' -@ {}'.format(max(1, len(self.pats) // 16))
            view_flags.append(v)
        if not view_flags:
            view_flags = None
        self.fast_merge_pats(view_flags)

    def compose_view_cmd(self, i, view_flags):
        if view_flags is None:
            view_cmd = ' <(gunzip -c'
        else:
            view_cmd = ' <({wt} cview {flags}'.format(wt=main_script, flags=view_flags[i])
            # view_cmd = ' <({wt} view {flags}'.format(wt=main_script, flags=view_flags[i])
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
                eprint(f'Invalid temp dir: {temp_dir}. Ignoring it')
            else:
                cmd += f' -T {temp_dir} '
        for i in range(len(self.pats)):
            cmd += self.compose_view_cmd(i, view_flags)
        cmd += ' | {} - '.format(collapse_pat_script)
        # cmd += f' | bedtools groupby -g 1-3 -c 4'  # TODO: check how many columns in self.pats[0] and use groupby instead of collapse_pat_script
        cmd += f' | bgzip > {self.outpath}'
        cmd = f'/bin/bash -c "{cmd}"'
        if self.args.verbose:
            eprint(cmd)
        subprocess.check_call(cmd, shell=True)

        if not op.isfile(self.outpath):
            raise IllegalArgumentError(f'[wt merge] Error: failed to create file {self.outpath}')
        Indxer(self.outpath).run()


def merge_betas(betas, opath, lbeta=False):
    """
    Merge all betas by summing their values element-wise, while keeping the dimensions
    :param betas: list of beta files
    :param opath: merged beta file
    """
    validate_file_list(betas)
    data = load_beta_data(betas[0]).astype(np.int)
    for b in betas[1:]:
        data += load_beta_data(b)

    # Trim / normalize to range [0, 256)
    data = trim_to_uint8(data, lbeta=lbeta)
    # Dump
    if lbeta and opath.endswith('.beta'):
        opath = opath[:-5] + '.lbeta'
    data.tofile(opath)
    return data



def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+')
    parser.add_argument('-p', '--prefix', help='Prefix of output file', required=True)
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if existed')
    parser.add_argument('-T', '--temp_dir', help='passed to "sort -m". Useful for merging very large pat files')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-l', '--lbeta', action='store_true', help='Use lbeta file (uint16) instead of beta (uint8)')
    parser.add_argument('--labels', nargs='+', help='labels for the mixed reads. '
                                                    'Default is None')
    add_GR_args(parser, bed_file=True)
    parser.add_argument('--min_len', type=int, default=None,
                        help='consider only reads covering at least MIN_LEN CpG sites [1]')
    parser.add_argument('--strict', action='store_true', help='Truncate reads that start/end outside the given region. '
                                                              'Only relevant if "region", "sites" '
                                                              'or "bed_file" flags are given.')

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

    # construct output path
    out_path = args.prefix + splitextgz(args.input_files[0])[1]

    if op.realpath(out_path) in [op.realpath(p) for p in args.input_files]:
        eprint('[wt merge] Error output path is identical ' \
                'to one of the input files {out_path}')
        return

    if not delete_or_skip(out_path, args.force):
        return

    files_type = splitextgz(input_files[0])[1][1:]

    if files_type in ('beta', 'bin'):
        merge_betas(input_files, out_path, args.lbeta)
    elif files_type == 'pat.gz':
        MergePats(input_files, args.prefix + '.pat.gz', args.labels, args).merge_pats()
    else:
        print('Unknown input format:', input_files[0])
        return


if __name__ == '__main__':
    main()
