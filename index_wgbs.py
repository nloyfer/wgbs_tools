#!/usr/bin/python3 -u

import argparse
import os.path as op
import sys
import subprocess
from utils_wgbs import IllegalArgumentError # todo: use this exception instead of returning non zero codes

def index_bgzipped_file(input_file):
    """
    Create a csi index for a pat file, assuming it's bgzipped
    :return: 0 iff succeeded
    """
    if input_file.endswith('unq.gz') or input_file.endswith('pat.gz'):
        tabix_params = ' -C -b 2 -e 2 '
    elif input_file.endswith('tsv.gz') or input_file.endswith('bed.gz'):
        tabix_params = ' -p bed '
    else:
        raise IllegalArgumentError('Couldn\'t index file of unknown format: {}'.format(input_file))

    cmd = 'tabix {} {}'.format(tabix_params, input_file)
    return subprocess.call(cmd, shell=True)


def bgzip_and_index_file(input_file):
    """
    bgzip and index uncompressed pat / unq input file
    :return: 0 iff succeeded in both operations
    """
    r = subprocess.call('bgzip -@ 4 ' + input_file, shell=True)
    if not r:
        r = index_bgzipped_file(input_file + '.gz')
    return r


def validate_file_suffix(input_file):
    """
    Make sure file suffix is one of the following:
    pat, pat.gz, unq, unq.gz
    :return: False it is not.
    """
    legit_file_suff = False
    for suf in ('.pat', '.pat.gz', '.unq', '.unq.gz'):
        if input_file.endswith(suf):
            legit_file_suff = True
    if not legit_file_suff:
        print('Index only supports pat or unq files.', input_file, file=sys.stderr)
        return False
    return True


def index_single_file(input_file, force):
    from utils_wgbs import delete_or_skip

    # validate input file
    if not (op.isfile(input_file)):
        print("no such file:", input_file, file=sys.stderr)
        return 1

    # # validate file suffix # todo: throws exception for tsv when loading blocks file in beta_vis.py
    # if not validate_file_suffix(input_file):
    #     return 1

    # if index already exists delete it or skip it
    if not delete_or_skip(input_file + '.csi', force):
        return 1

    if input_file.endswith('.gz'):

        # try indexing it:
        r = index_bgzipped_file(input_file)

        if r == 1:  # couldn't index because the file is gzipped instead of bgzipped
            # ungzip and bgzip
            subprocess.call(['unpigz', input_file])
            r = bgzip_and_index_file(op.splitext(input_file)[0])
            if not r:
                print('Success in indexing', file=sys.stderr)

    else:  # uncompressed file
        r = bgzip_and_index_file(input_file)

    if r:
        print('Failed with error code', r, file=sys.stderr)
    return r


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+', help='One or more file with extensions .pat[.gz] or .unq[.gz]')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing index file (csi) if existed')
    return parser.parse_args()


def main():
    """
    bgzip and index pat or unq files.
    Accepts single or multiple files.
    Files may be in various states: bgzipped, gzipped or uncompressed
    (i.e extensions .pat[.gz] or .unq[.gz]).
    bgzips them and generate an index for them (csi).
    """
    args = parse_args()

    for input_file in args.input_files:
        index_single_file(input_file, args.force)
    return 0


if __name__ == '__main__':
    main()

