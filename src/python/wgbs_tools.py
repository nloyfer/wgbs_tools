#!/usr/bin/env python3

import sys
import argparse
import importlib
from unittest.mock import patch

commands = [
    # view data
    'vis',
    'view',
    'cview',
    'convert',

    # convert beta to other formats
    'beta_to_blocks',
    'beta_to_table',
    'beta2bed',
    'beta2bw',
    'beta_cov',
    'beta_to_450k',

    # generate pats and betas
    'init_genome',
    'set_default_ref',
    'bam2pat',
    'index',
    'pat2beta',
    'bed2beta',
    'mix_pat',
    'merge',

    'segment',

    # less important
    'compare_betas',
    'homog',
    'find_markers',
    'bam_cpg_counts',
    'frag_len',
    'allele_split',
    'test_bimodal'
]

def main():
    if len(sys.argv) < 2 or (len(sys.argv) == 2 and sys.argv[1] in ('-h', '--help')):
        print_help()
        return

    parser = argparse.ArgumentParser(
        description='WGBS data analysis',
        usage='wgbstools <command> [<args>]')
    parser.add_argument('command', help='Subcommand to run')
    args = parser.parse_args(sys.argv[1:2])
    try:
        with patch.object(sys, 'argv', sys.argv[1:]):
            importlib.import_module(args.command).main()

    except ModuleNotFoundError as e:
        print_invalid_command(args.command)
        print_help()

    except ValueError as e:
        eprint(f'Invalid input argument\n{e}')
        return 1


# todo:
# bam2pat: Add reports to log file / stderr. e.g: % success, # sites covered, # reads extracted etc.
# pip install?

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def print_invalid_command(command):
    eprint('Invalid command:', f'\033[01;31m{command}\033[00m')
    from difflib import get_close_matches
    closets = [x for x in get_close_matches(command, commands)]
    if closets:
        eprint(f'did you mean \033[01;32m{closets[0]}\033[00m?')

def print_help(command=None):
    msg = '\nUsage: wgbstools <command> [<args>]'
    msg += '\nrun wgbstools <command> -h for more information'
    msg += '\nOptional commands:\n'
    print(msg)
    print(*commands, sep='\n')
    return 1


if __name__ == '__main__':
    main()

