#!/usr/bin/python3 -u

from unittest.mock import patch
import sys
from utils_wgbs import IllegalArgumentError
from view import main as view_main
from merge import main as merge_main
from index_wgbs import main as index_main
from vis import main as vis_main
from pat2beta import main as pat2beta_main
from cmp_betas import main as compare_beta_main
from beta_cov import main as beta_cov_main
from bed2beta import main as bed2beta_main
from beta2bed import main as beta2bed_main
from beta2bw import main as beta2bw_main
from beta_to_450k import main as beta_to_450k_main
from beta_to_blocks import main as beta_to_blocks_main
from bam2pat import main as bam2pat_main
from mix_pat import main as mix_pat_main
from dmb import main as dmb_main
import time
from datetime import timedelta


"""
Dependencies:
samtools
awk

python3.5
numpy
pandas
"""


callbacks = {
    'vis': vis_main,
    'view': view_main,
    'merge': merge_main,
    'index': index_main,
    'pat2beta': pat2beta_main,
    'compare_betas': compare_beta_main,
    'beta_cov': beta_cov_main,
    'bed2beta': bed2beta_main,
    'beta2bed': beta2bed_main,
    'beta2bw': beta2bw_main,
    'beta_to_450k': beta_to_450k_main,
    'beta_to_blocks': beta_to_blocks_main,
    'bam2pat': bam2pat_main,
    'mix_pat': mix_pat_main,
    'dmb': dmb_main
    # todo: unq2beta, collapse_beta_to_blocks
}

# todo:
# tests
# reference setup
# merge unq
# view -s unq (subsample unq)
# Default output directories: change from '.' to input src dir?
# Add reports to log file / stderr. e.g: % success, # sites covered, # reads extracted etc.
# translate region <-> sites - print nicely, with commas optional
# change pat to unq, and forget about olt pat format. show it with view?
# Change wgbs_tools.py to new name, update the print_help method.


def print_help(short=False):
    msg = 'Usage: wgbs_tools.py COMMAND [OPTIONS]\n\nOptional commands:\n'
    for key in sorted(callbacks.keys()):
        docs = callbacks[key].__doc__
        msg += '\n- ' + key
        if docs and not short:
            msg += docs
    if short:
        msg += '\nUse [-h] or COMMAND -h flag for additional information'
    print(msg, file=sys.stderr)


def main():

    if len(sys.argv) < 2 or (len(sys.argv) == 2 and sys.argv[1] in ('-h', '--help')):
        print_help()
        return

    print_time = False
    if '--time' in sys.argv:
        sys.argv.remove('--time')
        print_time = True
        start_time = time.time()
    try:
        command = sys.argv[1]
        if command not in callbacks.keys():
            print('Invalid command:', command, file=sys.stderr)
            print_help(short=True)
            return 1

        with patch.object(sys, 'argv', sys.argv[1:]):
            callbacks[command]()

    except IllegalArgumentError as e:
        print('Invalid input argument\n{}'.format(e), file=sys.stderr)
        return 1

    if print_time:
        print('time:', timedelta(seconds=time.time() - start_time), file=sys.stderr)


if __name__ == '__main__':
    main()
