#!/usr/bin/python3 -u

import os
import sys
import os.path as op
import argparse
import subprocess
from utils_wgbs import IllegalArgumentError, DIR


def eprint(*args, **kargs):
    print(*args, file=sys.stderr, **kargs)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', '-v', action='store_true')
    return parser.parse_args()


def compile(cmd, name, verbose):
    eprint('Compiling {}...'.format(name))
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.returncode:
        eprint("Failed compilation.\nCommand: {}\nreturn code: {}\nstderr:\n{}\nstdout:\n{}".
               format(cmd, p.returncode, output.decode(), error.decode()))
        raise RuntimeError('Failed compiling {}'.format(name))
    elif verbose:
        eprint(cmd)
        eprint(output.decode())
    eprint('success')

def main():
    args = parse_args()
    curdir = os.getcwd()
    try:
        os.chdir(DIR)
        compile_all(args)
    except RuntimeError as e:
        eprint(e)
        os.chdir(curdir)


def compile_all(args):
    # compile C++ files

    # stdin2beta (pat2beta)
    cmd = 'g++ -std=c++11 src/pat2beta/stdin2beta.cpp -o src/pat2beta/stdin2beta'
    compile(cmd, 'stdin2beta', args.verbose)

    # pat_sampler (pat view)
    cmd = 'g++ -std=c++11 src/pat_sampler/sampler.cpp -o src/pat_sampler/pat_sampler'
    compile(cmd, 'pat_sampler', args.verbose)

    # patter (bam2pat)
    cmd = 'g++ -std=c++11 pipeline_wgbs/patter.cpp -o pipeline_wgbs/patter'
    compile(cmd, 'patter', args.verbose)

    # match_maker (bam2pat)
    cmd = 'g++ -std=c++11 pipeline_wgbs/match_maker.cpp -o pipeline_wgbs/match_maker'
    compile(cmd, 'match_maker', args.verbose)


if __name__ == '__main__':
    main()
