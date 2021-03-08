#!/usr/bin/python3 -u
#!/usr/bin/env python3


import os
import sys
import matplotlib
if 'DISPLAY' not in os.environ.keys():
    matplotlib.use('Agg')
from collections import OrderedDict
from unittest.mock import patch
from utils_wgbs import IllegalArgumentError, eprint
from view import main as view_main
from cview import main as cview_main
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
from bamAddMethylData import main as bamAddMethylData_main
from mix_pat import main as mix_pat_main
from init_genome_ref_wgbs import main as init_genome_main
from frag_len import main as frag_len_main
from convert import main as convert_main
from segment_betas import main as segment_beta_main

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
    ('beta2bed', beta2bed_main),
    ('beta2bw', beta2bw_main),
    ('beta_cov', beta_cov_main),
    ('beta_to_450k', beta_to_450k_main),

    # generate pats and betas
    ('init_genome', init_genome_main),
    ('bam2pat', bam2pat_main),
    ('index', index_main),
    ('pat2beta', pat2beta_main),
    ('bed2beta', bed2beta_main),
    ('mix_pat', mix_pat_main),
    ('merge', merge_main),

    ('segment', segment_beta_main),

    # less important
    ('compare_betas', compare_beta_main),
    ('addMD2Bam', bamAddMethylData_main),
    ('frag_len', frag_len_main),
]
callbacks = OrderedDict(callbacks)

# todo:
# tests
# bam2pat: Add reports to log file / stderr. e.g: % success, # sites covered, # reads extracted etc.
# convert: translate region <-> sites - optional parsable format  
# Change wgbs_tools.py to new name, update the print_help method.


def print_help(short=False):
    msg = 'Usage: wgbs_tools.py COMMAND [OPTIONS]\n\nOptional commands:\n'
    for key in callbacks.keys():
        docs = callbacks[key].__doc__
        msg += '\n- ' + key
        if docs and not short:
            msg += docs
    if short:
        msg += '\nUse [-h] or COMMAND -h flag for additional information'
    eprint(msg)


def main():
    if len(sys.argv) < 2 or (len(sys.argv) == 2 and sys.argv[1] in ('-h', '--help')):
        print_help()
        return

    try:
        command = sys.argv[1]
        if command not in callbacks.keys():
            eprint('Invalid command:', command)
            print_help(short=True)
            return 1

        with patch.object(sys, 'argv', sys.argv[1:]):
            callbacks[command]()

    except IllegalArgumentError as e:
        eprint('Invalid input argument\n{}'.format(e))
        return 1


if __name__ == '__main__':
    main()
