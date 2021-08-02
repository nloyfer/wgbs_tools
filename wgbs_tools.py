#!/usr/bin/python3 -u
#!/usr/bin/env python3


import os
import sys
import matplotlib
if 'DISPLAY' not in os.environ.keys():
    matplotlib.use('Agg')
from collections import OrderedDict
from unittest.mock import patch
from src.python.utils_wgbs import IllegalArgumentError, eprint
from src.python.view import main as view_main
from src.python.cview import main as cview_main
from src.python.merge import main as merge_main
from src.python.index_wgbs import main as index_main
from src.python.vis import main as vis_main
from src.python.pat2beta import main as pat2beta_main
# from pat2pairs import main as pat2pairs_main
from src.python.cmp_betas import main as compare_beta_main
from src.python.beta_cov import main as beta_cov_main
from src.python.bed2beta import main as bed2beta_main
from src.python.beta2bed import main as beta2bed_main
from src.python.beta2bw import main as beta2bw_main
from src.python.beta_to_450k import main as beta_to_450k_main
from src.python.beta_to_blocks import main as beta_to_blocks_main
from src.python.betas_to_table import main as betas_to_table_main
from src.python.bam2pat import main as bam2pat_main
from src.python.bamAddMethylData import main as bamAddMethylData_main
from src.python.mix_pat import main as mix_pat_main
from src.python.init_genome_ref_wgbs import main as init_genome_main
from src.python.frag_len import main as frag_len_main
from src.python.convert import main as convert_main
from src.python.segment_betas import main as segment_beta_main
from src.python.homog import main as homog_main
# from dmb import main as dmb_main

"""
Dependencies:
samtools
awk

python3
numpy
pandas
"""

callbacks = [
    # view data
    ('vis', vis_main),
    ('view', view_main),
    ('cview', cview_main),
    ('convert', convert_main),

    # convert beta to other formats
    ('beta_to_blocks', beta_to_blocks_main),
    ('beta_to_table', betas_to_table_main),
    ('beta2bed', beta2bed_main),
    ('beta2bw', beta2bw_main),
    ('beta_cov', beta_cov_main),
    ('beta_to_450k', beta_to_450k_main),

    # generate pats and betas
    ('init_genome', init_genome_main),
    ('bam2pat', bam2pat_main),
    ('index', index_main),
    ('pat2beta', pat2beta_main),
    # ('pat2pairs', pat2pairs_main),
    ('bed2beta', bed2beta_main),
    ('mix_pat', mix_pat_main),
    ('merge', merge_main),

    ('segment', segment_beta_main),

    # less important
    ('compare_betas', compare_beta_main),
    ('homog', homog_main),
    # ('dmb', dmb_main),
    ('bam_cpg_counts', bamAddMethylData_main),
    ('frag_len', frag_len_main),
]
callbacks = OrderedDict(callbacks)

# todo:
# tests
# auto download genome (hg19, mm9)
# move all of the *py files to src directory
# pip install?
# bam2pat: Add reports to log file / stderr. e.g: % success, # sites covered, # reads extracted etc.
# convert: translate region <-> sites - optional parsable format  
# Change wgbs_tools.py to new name, update the print_help method.
# beta_to_450K: use formal annotations (now refering to cbio)


def print_help(short=False):
    msg = '\nUsage: wgbs_tools.py COMMAND [OPTIONS]\n\nOptional commands:\n'
    for key in callbacks.keys():
        docs = callbacks[key].__doc__
        msg += '\n- ' + key
        if docs and not short:
            msg += docs
    if short:
        msg += '\nUse [-h] or COMMAND -h flag for more information'
    eprint(msg)
    return 1


def print_closest(command):
    eprint('Close possible commands:')
    from difflib import get_close_matches
    for c in get_close_matches(command, callbacks.keys()):
        r = f'\033[01;32m{c}\033[00m'
        eprint(r)


def main():
    if len(sys.argv) < 2 or (len(sys.argv) == 2 and sys.argv[1] in ('-h', '--help')):
        return print_help()

    try:
        command = sys.argv[1]
        if command not in callbacks.keys():
            eprint('Invalid command:', f'\033[01;31m{command}\033[00m')
            print_closest(command)
            print_help(short=True)
            return 1

        with patch.object(sys, 'argv', sys.argv[1:]):
            callbacks[command]()

    except IllegalArgumentError as e:
        eprint(f'Invalid input argument\n{e}')
        return 1


if __name__ == '__main__':
    main()
