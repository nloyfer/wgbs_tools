#!/usr/bin/python3 -u

import argparse
import os.path as op
import sys
import subprocess
import multiprocessing
from utils_wgbs import delete_or_skip, splitextgz, IllegalArgumentError, eprint, add_multi_thread_args

# TODO: use inheritence
# TODO: drop tsv
class Pat:
    def __init__(self, input_file):
        self.input_file = input_file
        self.suffixes = ['pat']
        self.tabix_flags = ' -C -b 2 -e 2 -m 12 '
        self.sort_flags = '-k2,2n'
        self.ind_suff = '.csi'
        self.tosort = True


class BedTsv:
    def __init__(self, input_file):
        self.input_file = input_file
        self.suffixes = ['bed', 'tsv']
        self.tabix_flags = ' -p bed '
        self.sort_flags = '-k1,1 -k2,2n'
        self.ind_suff = '.tbi'
        self.tosort = False


class Indxer:
    def __init__(self, input_file, force=True, threads=multiprocessing.cpu_count()):
        self.force = force
        self.threads = threads
        self.in_file = input_file
        self.suff = splitextgz(self.in_file)[1][1:]
        c = BedTsv if 'bed' in self.suff or 'tsv' in self.suff else Pat
        self.ftype = c(input_file)
        self.validate_file()

    def index_bgzipped_file(self):
        """
        Create a csi index for the file, assuming it's bgzipped
        :return: 0 iff succeeded
        """
        cmd = 'tabix {} {}.gz'.format(self.ftype.tabix_flags, self.in_file)
        r = subprocess.check_call(cmd, shell=True, stderr=subprocess.PIPE)
        return r

    def bgzip(self):
        """ bgzip uncompressed input file """
        f = self.in_file
        # check if file is sorted. If not, sort it:
        flags = self.ftype.sort_flags
        if self.ftype.tosort and 'pat' in self.suff:
            if subprocess.call(f'sort {flags} -cs {f}', shell=True):
                eprint(f'{f} is not sorted. Sorting...')
                subprocess.check_call('sort {fl} {f} -o {f}'.format(fl=flags, f=f), shell=True)

        # bgzip the file:
        cmd = f'bgzip -@ {self.threads} -f {f}'
        subprocess.check_call(cmd, shell=True)

    def validate_file(self):
        """
        Make sure file exists, and its suffix is one of the following:
        pat, bed, tsv [.gz]
        """

        if not op.isfile(self.in_file):
            raise IllegalArgumentError("no such file: {self.in_file}")

        suffs = self.ftype.suffixes
        if not self.suff in [x + '.gz' for x in suffs] + suffs:
            raise IllegalArgumentError('Index only supports pat, bed, tsv formats')

    def run(self):

        # if index already exists delete it or skip it
        if not delete_or_skip(self.in_file + self.ftype.ind_suff, self.force):
            return

        if self.in_file.endswith('.gz'):
            self.in_file = op.splitext(self.in_file)[0]
            # try indexing it:
            if not self.index_bgzipped_file():
                return  # success
            # couldn't index because the file is gzipped instead of bgzipped
            subprocess.check_call(['gunzip', self.in_file + '.gz'])

        self.bgzip()
        self.index_bgzipped_file()
#        print('success')


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+', help='One or more file with extensions .pat[.gz] or .bed[.gz]')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing index file (csi) if existed')
    add_multi_thread_args(parser)
    return parser.parse_args()


def main():
    """
    bgzip and index pat or bed files.
    Accepts single or multiple files.
    Files may be in various states: bgzipped, gzipped or uncompressed
    (i.e extensions .pat[.gz] or .bed[.gz]).
    bgzips them and generate an index for them (csi).
    """
    args = parse_args()
    for input_file in args.input_files:
        Indxer(input_file, args.force, args.threads).run()
    return 0


if __name__ == '__main__':
    main()

