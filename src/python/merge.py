#!/usr/bin/python3 -u

import argparse
import subprocess
import os.path as op
from pathlib import Path
import numpy as np
from index import Indxer
from utils_wgbs import validate_file_list, splitextgz, delete_or_skip, \
        trim_to_uint8, load_beta_data, collapse_pat_script, \
        IllegalArgumentError, main_script, eprint, pretty_name, \
        GenomeRefPaths, safe_remove
from cview import add_view_flags
from genomic_region import GenomicRegion


def validate_labels(labels, pats, required=False):
    if not labels and not required:
        return labels
    if not labels:
        labels = [pretty_name(p) for p in pats]

    if len(labels) != len(pats):
        raise IllegalArgumentError('[wt mix] len(labels) != len(files)')
    if len(labels) != len(set(labels)):
        raise IllegalArgumentError('[wt mix] duplicated labels')
    return labels


def extract_view_flags(args):
    v = f' --genome {args.genome}'
    if args.strict:
        v += ' --strict'
    if args.strip:
        v += ' --strip'
    if args.min_len:
        v += f' --min_len {args.min_len}'
    if args.bed_file is not None:
        v += f' -L {args.bed_file}'
    gr = GenomicRegion(args)
    if not gr.is_whole():
        v += ' -s {}-{}'.format(*gr.sites)
    return v


class MergePats:
    def __init__(self, pats, outpath, labels, args):
        self.args = args
        self.gr = GenomicRegion(args)
        self.pats = pats
        validate_file_list(self.pats, force_suff='.pat.gz')
        self.outpath = outpath
        self.labels = validate_labels(labels, pats)

    def merge_pats(self):
        view_flags = [extract_view_flags(self.args)] * len(self.pats)

        if self.gr.is_whole():
            self.merge_by_chrom(view_flags)
        else:
            self.fast_merge_pats(view_flags)
        Indxer(self.outpath).run()

    def compose_view_cmd(self, i, view_flags):
        if not view_flags:
            view_cmd = ' <(gunzip -c'
        else:
            view_cmd = f' <({main_script} cview {view_flags[i]}'
        view_cmd += f' {self.pats[i]}'
        tagcmd = ''
        if self.labels is not None:
            tagcmd = f' | sed s/\$/\'\t\'{self.labels[i]}/'
        view_cmd += f'{tagcmd})'
        return view_cmd

    def fast_merge_pats(self, view_flags=None, append=False):
        """ Use piping and sort -m to merge pat files w/o intermediate files """
        cmd = 'sort -m -k2,2n -k3,3'

        # add/validate temp dir
        if self.args.temp_dir:
            temp_dir = self.args.temp_dir
            if not op.isdir(temp_dir):
                eprint(f'Invalid temp dir: {temp_dir}. Ignoring it')
            else:
                cmd += f' -T {temp_dir} '

        # compose view command for each pat
        for i in range(len(self.pats)):
            cmd += self.compose_view_cmd(i, view_flags)
        cmd += f' | {collapse_pat_script} - '
        # cmd += f' | bedtools groupby -g 1-3 -c 4'  # TODO: check how many columns in self.pats[0] and use groupby instead of collapse_pat_script
        if append:
            cmd += f' | bgzip >> {self.outpath}'
        else:
            cmd += f' | bgzip > {self.outpath}'
        cmd = f'/bin/bash -c "{cmd}"'
        if self.args.verbose:
            eprint(cmd)
        subprocess.check_call(cmd, shell=True)

        if not op.isfile(self.outpath):
            raise IllegalArgumentError(f'[wt merge] Error: failed to create file {self.outpath}')

    def merge_by_chrom(self, view_flags):
        # merge pat files chrom by chrom when we merge whole-genome (much faster because it saves the sorting process lots of time)
        sorted_chroms = GenomeRefPaths(self.args.genome).get_chroms()

        # dump merged pat chromosome by chromosome,
        # so that the unix sort will run on each chromosome separately
        safe_remove(self.outpath)
        Path(self.outpath).touch()
        if self.args.verbose:
            eprint('[wt merge] Merging chrom by chrom...')
        for chrom in sorted_chroms:
            if self.args.verbose:
                eprint(f'[wt merge] {chrom}')
            chrom_view_flags = [v + f' -r {chrom}' for v in view_flags]
            self.fast_merge_pats(chrom_view_flags, append=True)


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
    parser = add_view_flags(parser, sub_sample=False, out_path=False)
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
                f'to one of the input files {out_path}')
        return

    if not delete_or_skip(out_path, args.force):
        return

    files_type = splitextgz(input_files[0])[1][1:]

    if files_type in ('beta', 'lbeta', 'bin'):
        merge_betas(input_files, out_path, args.lbeta)
    elif files_type == 'pat.gz':
        MergePats(input_files, args.prefix + '.pat.gz', args.labels, args).merge_pats()
    else:
        print('Unknown input format:', input_files[0])
        return


if __name__ == '__main__':
    main()
