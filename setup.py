#!/usr/bin/python3 -u

import os
import sys
import os.path as op
import argparse
import subprocess
from utils_wgbs import DIR, eprint, validate_single_file, IllegalArgumentError


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', '-v', action='store_true')
    return parser.parse_args()


def compile_single(cmd, name, verbose):
    """ Compile a single c++ module """
    eprint(f'Compiling {name}...')

    # subprocess g++ command
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()

    # If failed:
    if p.returncode:
        eprint("Failed compilation.\nCommand: {}\nreturn code: {}\nstderr:\n{}\nstdout:\n{}".
               format(cmd, p.returncode, output.decode(), error.decode()))
        eprint(f'Failed compiling {name}')
        # raise RuntimeError('Failed compiling {}'.format(name))

    # If succeeded:
    elif verbose:
        eprint(cmd)
        eprint(output.decode())
    eprint('success')


def compile_all(args):
    curdir = os.getcwd()
    try:
        os.chdir(DIR)

        # compile C++ files

        # stdin2beta (pat2beta)
        cmd = 'g++ -std=c++11 src/pat2beta/stdin2beta.cpp -o src/pat2beta/stdin2beta'
        compile_single(cmd, 'stdin2beta', args.verbose)

        # pat_sampler (pat view)
        cmd = 'g++ -std=c++11 src/pat_sampler/sampler.cpp -o src/pat_sampler/pat_sampler'
        compile_single(cmd, 'pat_sampler', args.verbose)

        # patter (bam2pat)
        cmd = 'g++ -std=c++11 src/pipeline_wgbs/patter.cpp -o src/pipeline_wgbs/patter'
        compile_single(cmd, 'patter', args.verbose)

        # match_maker (bam2pat)
        cmd = 'g++ -std=c++11 src/pipeline_wgbs/match_maker.cpp -o src/pipeline_wgbs/match_maker'
        compile_single(cmd, 'match_maker', args.verbose)

        # segmentor (segment)
        cmd = 'g++ -std=c++11 src/segment_betas/main.cpp src/segment_betas/segmentor.cpp -o src/segment_betas/segmentor '
        compile_single(cmd, 'segmentor', args.verbose)

        # cview (cview)
        cmd = 'g++ -c -o src/cview/main.o src/cview/main.cpp -std=c++11; g++ -c -o src/cview/cview.o src/cview/cview.cpp ; g++ -o src/cview/cview src/cview/main.o src/cview/cview.o -std=c++11'
        compile_single(cmd, 'cview', args.verbose)

    except RuntimeError as e:
        eprint(e)
        os.chdir(curdir)
        return

    os.chdir(curdir)


# def validate_index_file(bedfile):
    # if not op.isfile(bedfile + '.tbi'):
        # eprint(f'No tbi found for file {bedfile}. Attempting to index it...')
        # from index_wgbs import Indxer
        # Indxer(bedfile).run()


# def validate_file(file):
    # if file is None:
        # return
    # validate_single_file(file)


# def config_file(args):
    # from config_wgbs import default_anno_file, default_blocks_path, ilmn2cpg_dict

    # validate_file(default_anno_file)
    # validate_index_file(default_anno_file)

    # validate_file(default_blocks_path)
    # validate_index_file(default_blocks_path)

    # validate_file(ilmn2cpg_dict)
    # return

# TODO:
# 1. blocks file: make a symbolic link to references/genome/blocks.bed.gz[.tbi]
#    Maybe even upload an ultra compressed (diffs) version of the latest hg19 blocks to GitHub.
# 2. add shebang line to wgbs_tools.py. Maybe add it to PATH, or add instructions in the README.md
# 3. Annotaions.gz file - same as blocks.

def main():
    args = parse_args()
    compile_all(args)
    # config_file(args)


if __name__ == '__main__':
    main()
