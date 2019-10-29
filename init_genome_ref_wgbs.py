#!/usr/bin/python3 -u

import argparse
from utils_wgbs import  delete_or_skip, validate_single_file, eprint
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
    print('Dumping data frame to path:', path)
    df.to_csv(path, index=None, header=None, sep='\t')


def bgzip_tabix_dict(dict_path):
    eprint('bgzip and index...')
    subprocess.check_call('bgzip -@ 4 -f ' + dict_path, shell=True)
    subprocess.check_call('tabix -Cf -b 2 -e 2 {}.gz'.format(dict_path), shell=True)


def init_ref_files(ref_fasta, out_dir):
    df = pd.DataFrame()     # Full dict: Chr    Start   End    CpGIndex
    c_cpg_sizes = {}        # chromosome sizes (#sites)
    c_sizes = {}            # chromosomes sizes (#bp)

    eprint('Loading fasta...')
    fiter = fasta_iter(ref_fasta)
    # Parse genome reference fasta, chromosome by chromosome
    for ff in sorted(fiter, key=pair_chromosome_order):
        chrName, seq = ff
        if not is_valid_chrome(chrName):
            # print('Invalid:', chrName)
            continue
        print(chrName)
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


def set_output_dir(name):
    out_dir = op.join(op.dirname(__file__), 'references')
    if not op.isdir(out_dir):
        os.mkdir(out_dir)
    return op.join(out_dir, name)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('genome_ref', help='path to a genome *.fa file')
    parser.add_argument('name', help='name of the genome (e.g. hg19, mm9...).\n'
                                     'A directory of this name will be created '
                                     'in references/')
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
    out_dir = set_output_dir(args.name)
    eprint('Setting up genome reference files in {}'.format(out_dir))
    if not delete_or_skip(op.join(out_dir, 'CpG.bed.gz'), args.force):
        return

    init_ref_files(genome_ref, out_dir)


if __name__ == '__main__':
    main()
