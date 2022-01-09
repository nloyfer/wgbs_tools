#!/usr/bin/python3 -u

import argparse
import os
import os.path as op
import sys
import subprocess as sp
import multiprocessing
import tempfile
from utils_wgbs import delete_or_skip, splitextgz, IllegalArgumentError, eprint, \
        add_multi_thread_args, validate_single_file, COORDS_COLS5

class Pat:
    def __init__(self):
        self.suff = 'pat'
        self.tabix_flags = ' -C -b 2 -e 2 -m 12 '
        self.sort_flags = '-k2,2n -k3,3'
        self.ind_suff = '.csi'


class Bed:
    def __init__(self):
        self.suff = 'bed'
        self.tabix_flags = ' -p bed '
        self.sort_flags = '-k4,4n'
        self.ind_suff = '.tbi'


def tabix_fai_workaround(in_file):
    # if tabix version is >1.9, We need to override their file type deduction
    # see https://github.com/samtools/htslib/issues/1347
    # Will probably be fixed in release 1.15 of htslib (TODO: update here when they release)

    try:
        # only relevant for files with 5 columns, where the 5'th is an integer
        with open(in_file, 'w') as f:
            tokens = f.readline().split('\t')
        if not (len(tokens) == 5 and tokens[4].isdigit()):
            return

        # only relevant for tabix version >1.9
        txt = sp.check_output('tabix --version', shell=True).decode().split()[2]
        try:
            tabix_version = float(txt)
        except ValueError as e:
            tabix_version = None
        if tabix_version is not None and tabix_version <= 1.9:
            return
    except:
        return

    # add a header line
    tmp_name = op.basename(in_file) + next(tempfile._get_candidate_names())
    temp_path = op.join(op.dirname(in_file), tmp_name)
    try:
        with open(temp_path, 'w') as f:
            f.write('#' + '\t'.join(COORDS_COLS5) + '\n')
        cmd = 'cat {i} >> {t} && mv -f {t} {i}'.format(i=in_file, t=temp_path)
        sp.check_call(cmd, shell=True)
    finally:
        if op.isfile(temp_path):
            os.remove(temp_path)

class Indxer:
    def __init__(self, input_file, force=True,
                 threads=multiprocessing.cpu_count()):
        self.force = force
        self.threads = threads
        self.in_file = input_file
        self.suff = splitextgz(self.in_file)[1][1:]
        self.ftype = Bed() if 'bed' in self.suff else Pat()
        self.validate_file()

    def is_file_gzipped(self):
        return self.in_file.endswith('.gz') and \
              'BGZF' not in sp.check_output(f'htsfile {self.in_file}',
                                            shell=True).decode()

    def index_bgzipped_file(self):
        """
        Create a csi index for the file, assuming it's bgzipped
        :return: 0 iff succeeded
        """
        cmd = 'tabix {} {}'.format(self.ftype.tabix_flags, self.in_file)
        return sp.check_call(cmd, shell=True, stderr=sp.PIPE)


    def bgzip(self):
        """ bgzip uncompressed input file """

        # check if file is sorted. If not, sort it:
        scmd = f'sort {self.ftype.sort_flags} {self.in_file}'
        if sp.call(scmd + ' -sc', shell=True):
            eprint(f'[wt index] {self.in_file} is not sorted. Sorting...')
            sp.check_call(scmd + f' -o {self.in_file}', shell=True)

        # workaround for some tabix versions
        tabix_fai_workaround(self.in_file)

        # bgzip the file:
        cmd = f'bgzip -@ {self.threads} -f {self.in_file}'
        sp.check_call(cmd, shell=True)
        self.in_file += '.gz'

    def validate_file(self):
        """
        Make sure file exists, and its suffix is one of the following:
        pat, bed [.gz]
        """

        if not op.isfile(self.in_file):
            raise IllegalArgumentError(f'no such file: {self.in_file}')

        suffs = self.ftype.suff
        if not self.suff in (suffs, suffs + '.gz'):
            raise IllegalArgumentError('Index only supports pat, bed formats')

    def run(self):

        # if index already exists delete it or skip it
        if not delete_or_skip(self.in_file + self.ftype.ind_suff, self.force):
            return

        # if file is gzipped instead of bgzipped, uncompress it
        if self.is_file_gzipped():
            sp.check_call(['gunzip', self.in_file])
            self.in_file = self.in_file[:-3]

        if not self.in_file.endswith('.gz'):
            self.bgzip()
        self.index_bgzipped_file()


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+',
            help=f'One or more file with extensions .pat[.gz] or .bed[.gz]. Bed files must have the columns {COORDS_COLS5}')
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


if __name__ == '__main__':
    main()

