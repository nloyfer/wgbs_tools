#!/usr/bin/env python3

import sys
import argparse
import importlib
from unittest.mock import patch

VERSION = '0.2.0'

commands = [
    # view data
    'vis',
    'view',
    'cview',
    'convert',
    'pat_fig',

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
    'add_cpg_counts',
    'frag_len',
    'split_by_allele',
    'split_by_meth',
    'test_bimodal'
]

def main():
    if len(sys.argv) < 2 or (len(sys.argv) == 2 and sys.argv[1] in ('-h', '--help')):
        print_help()
        return
    if '--version' in sys.argv:
        print('wgbstools version', VERSION)
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
        if args.command not in str(e):
            raise e
        print_invalid_command(args.command)
        print_help()

    except ValueError as e:
        eprint(f'Invalid input argument\n{e}')
        return 1


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def print_invalid_command(command):
    eprint('Invalid command:', f'\033[01;31m{command}\033[00m')
    from difflib import get_close_matches
    closets = get_close_matches(command, commands)
    if closets:
        eprint(f'did you mean \033[01;32m{closets[0]}\033[00m?')

def print_help():
    msg = '\nUsage: wgbstools <command> [<args>]'
    msg += '\nrun wgbstools <command> -h for more information'
    msg += '\nOptional commands:\n'
    print(msg)
    print(*commands, sep='\n')
    return 1


if __name__ == '__main__':
    main()
