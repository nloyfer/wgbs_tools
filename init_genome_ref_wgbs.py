#!/usr/bin/python3 -u

import argparse
from utils_wgbs import validate_files_list, splitextgz, delete_or_skip, trim_to_uint8, load_beta_data, \
    validate_single_file
import subprocess
import numpy as np
import os.path as op
import sys
import os
from index_wgbs import index_single_file
import re
import numpy as np
import pandas as pd
import os.path as op
import os
from itertools import groupby
import subprocess
import sys


gref = '/cs/cbio/netanel/indexes/WholeGenomeFasta/genome.fa'


def fasta_iter(fasta_name):
    """ Based on: https://www.biostars.org/p/710/ """
    faiter = (x[1] for x in groupby(open(fasta_name), lambda line: line[0] == ">"))
    for header in faiter:
        yield (header.__next__()[1:].strip(), "".join(s.strip() for s in faiter.__next__()))


def pair_chromosome_order(pair):
    return chromosome_order(pair[0])


def chromosome_order(c):
    if not c.startswith('chr'):
        raise RuntimeError('Invalid chromosome' + c)
    c = c[3:]
    if c.isdigit():
        return int(c)
    elif c == 'X':
        return 100
    elif c == 'Y':
        return 101
    elif c == 'M':
        return 102
    return 103


def is_valid_chrome(chrome):
    return re.match(r'^chr([\d]+|[XYM])$', chrome)


def dump_df(df, path):
    df.to_csv(path, index=None, header=None, sep='\t')


def bgzip_tabix_dict(dict_path):
    print('bgzip and index...', file=sys.stderr)
    r = subprocess.call('bgzip -@ 4 -f ' + dict_path, shell=True)
    if r:
        print('Error while bgzip', file=sys.stderr)
        return
    r = subprocess.call('tabix -Cf -b 2 -e 2 {}.gz'.format(dict_path), shell=True)
    if r:
        print('Error while indexing with tabix', file=sys.stderr)   # todo: make exception


def init_ref_files(ref_fasta, out_dir):
    df = pd.DataFrame()     # Full dict: Chr    Start   End    CpGIndex
    c_cpg_sizes = {}        # chromosome sizes (#sites)
    c_sizes = {}            # chromosomes sizes (#bp)

    print('Loading fasta...', file=sys.stderr)
    fiter = fasta_iter(ref_fasta)
    # Parse genome reference fasta, chromosome by chromosome
    for ff in sorted(fiter, key=pair_chromosome_order):
        chrName, seq = ff
        if not is_valid_chrome(chrName):
            # print('Invalid:', chrName)
            continue
        # else:
        #     print('valid:', chrName)
        tf = pd.DataFrame([m.start() + 1 for m in re.finditer('CG', seq.upper())], columns=['loc'])
        tf['chr'] = chrName
        df = pd.concat([df, tf[['chr', 'loc']]])

        # sizes
        c_sizes[chrName] = len(seq)
        c_cpg_sizes[chrName] = tf.shape[0]

    df.reset_index(drop=True, inplace=True)
    df['site'] = df.index + 1
    dict_path = op.join(out_dir, 'CpG.bed')
    dump_df(df, dict_path)

    # bgzip and tabix
    bgzip_tabix_dict(dict_path)

    # reverse dict
    np.array(df['loc'] + 1, dtype=np.uint32).tofile(dict_path.replace('.bed', '.rev.bin'))

    # chromosomes sizes
    c_sizes = pd.DataFrame.from_dict(c_sizes, orient='index').reset_index()
    c_cpg_sizes = pd.DataFrame.from_dict(c_cpg_sizes, orient='index').reset_index()

    c_size_path = op.join(out_dir, 'chrome.size')
    c_cpg_size_path = op.join(out_dir, 'CpG.chrome.size')

    dump_df(c_sizes, c_size_path)
    dump_df(c_cpg_sizes, c_cpg_size_path)


# def main1():
#     # outpath = 'a.tsv'
#     # gref = 'test.fa'
#     out_dir = 'references'
#     test(gref, out_dir)
#     # if gref != '/cs/cbio/netanel/indexes/WholeGenomeFasta/genome.fa':
#     #     os.system('cat ' + outpath)
#     # GenomeRefParser(gref, 30).read_ref()


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('genome_ref')
    parser.add_argument('-o', '--out_dir')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')
    args = parser.parse_args()
    return args


def main():
    """
    Init genome reference.
    """
    args = parse_args()

    # validate input files
    genome_ref = args.genome_ref
    validate_single_file(genome_ref, '.fa')

    # construct output_dir
    out_dir = args.out_dir
    if not out_dir:
        out_dir = op.dirname(genome_ref)
    if not op.isdir(out_dir):
        print('Invalid output dir:', out_dir, file=sys.stderr)
        return
    if not delete_or_skip(op.join(out_dir, 'CpG.bed.gz'), args.force):
        return

    init_ref_files(genome_ref, out_dir)

if __name__ == '__main__':
    main()
