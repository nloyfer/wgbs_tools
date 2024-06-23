#!/usr/bin/python3 -u

import os
import sys
import os.path as op
import argparse
import subprocess
from src.python.utils_wgbs import eprint, validate_single_file, IllegalArgumentError


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--targets', '-t', nargs='+', help='compile only these modules')
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
        eprint('\033[01;31mFAIL\033[00m')
        eprint("Failed compilation.\nCommand: {}\nreturn code: {}\nstderr:\n{}\nstdout:\n{}".
               format(cmd, p.returncode, output.decode(), error.decode()))
        eprint(f'Failed compiling {name}')
        # raise RuntimeError('Failed compiling {}'.format(name))

    # If succeeded:
    elif verbose:
        eprint(cmd)
        eprint(output.decode())
    else:
        eprint('\033[01;32mSUCCESS\033[00m')

modules = {
        'stdin2beta'    : 'g++ -std=c++11 src/pat2beta/stdin2beta.cpp -o src/pat2beta/stdin2beta',
        'mask_pat'      : 'g++ -std=c++11 -c -o src/pat2beta/mask_pat.o src/pat2beta/mask_pat.cpp && '\
                          'g++ -std=c++11 -c -o src/pipeline_wgbs/patter_utils.o src/pipeline_wgbs/patter_utils.cpp && '\
                          'g++ -std=c++11 -o src/pat2beta/mask_pat src/pat2beta/mask_pat.o src/pipeline_wgbs/patter_utils.o',
        # 'stdin2pairs'   : 'g++ -std=c++11 src/pat2beta/stdin2pairs.cpp -o src/pat2beta/stdin2pairs',
        'pat_sampler'   : 'g++ -std=c++11 src/pat_sampler/sampler.cpp -o src/pat_sampler/pat_sampler',
        'patter'        : 'g++ -std=c++11 -c -o src/pipeline_wgbs/main.o src/pipeline_wgbs/main.cpp && '\
                          'g++ -std=c++11 -c -o src/pipeline_wgbs/patter.o src/pipeline_wgbs/patter.cpp && '\
                          'g++ -std=c++11 -c -o src/pipeline_wgbs/patter_utils.o src/pipeline_wgbs/patter_utils.cpp && '\
                          'g++ -std=c++11 -c -o src/pipeline_wgbs/ont.o src/pipeline_wgbs/ont.cpp && '\
                          'g++ -std=c++11 -o src/pipeline_wgbs/patter src/pipeline_wgbs/main.o src/pipeline_wgbs/patter_utils.o src/pipeline_wgbs/patter.o src/pipeline_wgbs/ont.o',
        'bp_patter'     : 'g++ -std=c++11 src/pipeline_wgbs/blueprint/patter.cpp -o src/pipeline_wgbs/blueprint/patter',
        'snp_patter'    : 'g++ -std=c++11 -c -o src/pipeline_wgbs/snp_patter.o src/pipeline_wgbs/snp_patter.cpp && '\
                          'g++ -std=c++11 -c -o src/pipeline_wgbs/patter_utils.o src/pipeline_wgbs/patter_utils.cpp && '\
                          'g++ -std=c++11 -o src/pipeline_wgbs/snp_patter src/pipeline_wgbs/snp_patter.o src/pipeline_wgbs/patter_utils.o',
        'match_maker'   : 'g++ -std=c++11 src/pipeline_wgbs/match_maker.cpp src/pipeline_wgbs/patter_utils.cpp -o src/pipeline_wgbs/match_maker',
        'segmentor'     : 'g++ -std=c++11 src/segment_betas/main.cpp src/segment_betas/segmentor.cpp -o src/segment_betas/segmentor ',
        'cview'         : 'g++ -std=c++11 -c -o src/cview/cview.o src/cview/cview.cpp && '\
                          'g++ -std=c++11 -c -o src/pipeline_wgbs/patter_utils.o src/pipeline_wgbs/patter_utils.cpp && '\
                          'g++ -std=c++11 -o src/cview/cview src/cview/cview.o src/pipeline_wgbs/patter_utils.o',
        'homog'         : 'g++ -std=c++11 -c -o src/homog/homog.o src/homog/homog.cpp && '\
                          'g++ -std=c++11 -c -o src/pipeline_wgbs/patter_utils.o src/pipeline_wgbs/patter_utils.cpp && '\
                          'g++ -std=c++11 -o src/homog/homog src/homog/homog.o src/pipeline_wgbs/patter_utils.o',
        'add_cpg_counts': 'g++ -std=c++11 -c -o src/pipeline_wgbs/add_cpg_counts.o src/pipeline_wgbs/add_cpg_counts.cpp && '\
                          'g++ -std=c++11 -c -o src/pipeline_wgbs/patter_utils.o src/pipeline_wgbs/patter_utils.cpp && '\
                          'g++ -std=c++11 -o src/pipeline_wgbs/add_cpg_counts src/pipeline_wgbs/add_cpg_counts.o src/pipeline_wgbs/patter_utils.o',
        'add_loci'      : 'g++ -std=c++11 -pthread src/cpg2bed/add_loci.cpp src/cpg2bed/cpg_dict.cpp -o src/cpg2bed/add_loci'
        }

def compile_all(args):
    curdir = os.getcwd()
    try:
        os.chdir(op.dirname(op.realpath(__file__)))

        # compile C++ files
        for module in modules.keys():
            if args.targets and module not in args.targets:
                continue
            compile_single(modules[module], module, args.verbose)

    except RuntimeError as e:
        eprint(e)
        os.chdir(curdir)
        return

    os.chdir(curdir)


def main():
    args = parse_args()
    compile_all(args)
    # config_file(args)


if __name__ == '__main__':
    main()
