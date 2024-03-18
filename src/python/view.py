#!/usr/bin/python3 -u

import argparse
import os.path as op
import numpy as np
from utils_wgbs import load_beta_data, validate_single_file, \
    IllegalArgumentError, catch_BrokenPipeError, view_beta_script, \
    view_lbeta_script, beta_sanity_check
from genomic_region import GenomicRegion
from cview import cview, subprocess_wrap_sigpipe, add_view_flags


# def is_bed_disjoint(b):
    # if b.endswith('.gz'):
        # return      # fail quietly to warn user
    # cmd = f"""/bin/bash -c 'diff {b} <(bedtools intersect -a {b} -b {b} -wa)' > /dev/null """
    # if subprocess.call(cmd, shell=True):
        # eprint(f'[wt view] WARNING: bed file {b} regions are not disjoint.')



####################
#                  #
#  Loading beta    #
#                  #
####################

def view_other_bin(bin_path, args):
    # view bin files. Minimal support. Works very slow for whole genome.
    gr = GenomicRegion(args)
    data = load_beta_data(bin_path, gr.sites)
    np.savetxt('/dev/stdout', data, fmt='%s', delimiter='\t')


def bview_build_cmd(beta_path, gr, bed_path):
    beta_sanity_check(beta_path, gr.genome)
    # compose a shell command to output a beta file to stdout
    if beta_path.endswith('.beta'):
        vs = view_beta_script
    elif beta_path.endswith('.lbeta'):
        vs = view_lbeta_script
    else:
        raise IllegalArgumentError('Invalid input format for view beta:', beta_path)

    cmd = f'{vs} {gr.genome.revdict_path} {beta_path} '
    if not gr.is_whole():
        cmd += f' {gr.chrom} {gr.sites[0]} {gr.nr_sites}'
    if bed_path:
        validate_single_file(bed_path)
        cmd += f' | bedtools intersect -b {bed_path} -a stdin -wa '
    return cmd


def view_beta(beta_path, gr, opath, bed_path):
    """
    View beta file in given region/sites range/s
    :param beta_path: beta file path
    :param gr: a GenomicRegion object
    :param opath: output path (or stdout)
    :param bed_path: path to a bed file
    """

    cmd = bview_build_cmd(beta_path, gr, bed_path)

    if opath is not None:
        if opath.endswith('.gz'):
            cmd += ' | gzip -c '
        cmd += f' > {opath}'
    subprocess_wrap_sigpipe(cmd)

##########################
#                        #
#         Main           #
#                        #
##########################


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_file')
    parser = add_view_flags(parser)
    return parser


def main():
    """
    View the content of input file (pat/beta) as plain text.
    Possible filter by genomic region or sites range
    Output to stdout as default
    """
    parser = parse_args()
    args = parser.parse_args()

    if args.sub_sample is not None and not 1 >= args.sub_sample >= 0:
        parser.error('[wt view] sub-sampling rate must be within [0.0, 1.0]')

    # validate input file
    input_file = args.input_file
    validate_single_file(input_file)


    try:
        if op.splitext(input_file)[1] in ('.lbeta', '.beta'):
            gr = GenomicRegion(args)
            view_beta(input_file, gr, args.out_path, args.bed_file)
        elif op.splitext(input_file)[1] in ('.bin',):
            view_other_bin(input_file, args)
        elif input_file.endswith('.pat.gz'):
            cview(input_file, args)
        else:
            raise IllegalArgumentError('Unknown input format:', input_file)

    except BrokenPipeError:
        catch_BrokenPipeError()


if __name__ == '__main__':
    main()
