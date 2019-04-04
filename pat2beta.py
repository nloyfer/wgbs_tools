#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_single_file, PAT2BETA_TOOL, delete_or_skip, splitextgz, IllegalArgumentError
import subprocess
import os.path as op


def pat2beta(pat_path, out_dir, force=True):
    if pat_path.endswith('.pat.gz'):
        cmd = 'gunzip -cd'
    elif pat_path.endswith('.pat'):
        cmd = 'cat'
    else:
        raise IllegalArgumentError('Invalid pat suffix: {}'.format(pat_path))

    out_beta = op.join(out_dir, splitextgz(op.basename(pat_path))[0] + '.beta')
    if not delete_or_skip(out_beta, force):
        return
    cmd += ' {} | {} {}'.format(pat_path, PAT2BETA_TOOL, out_beta)
    r = subprocess.call(cmd, shell=True)
    if r:
        raise IllegalArgumentError('Failed generating beta from pat. error code: {}. {}'.format(r, pat_path))
    return out_beta


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat_path', help='A pat[.gz] file')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if existed')
    parser.add_argument('-o', '--out_dir', help='Output directory for the beta file. [.]', default='.')
    # todo: add genome paramter and support mm9
    return parser.parse_args()


def main():
    """
    Generate a beta file from a pat file
    """
    args = parse_args()
    validate_single_file(args.pat_path)
    pat2beta(args.pat_path, args.out_dir, args.force)


if __name__ == '__main__':
    main()

