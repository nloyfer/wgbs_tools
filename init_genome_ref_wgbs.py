#!/usr/bin/python3 -u

import argparse
from utils_wgbs import delete_or_skip, validate_single_file, eprint, IllegalArgumentError, DIR
import re
import numpy as np
import pandas as pd
import os.path as op
import os
from itertools import groupby
import subprocess
from multiprocessing import Pool
import multiprocessing


class InitGenome:
    def __init__(self, args):
        self.args = args
        self.ref_path = args.genome_ref
        self.force = args.force
        self.name = args.name
        self.out_dir = self.make_output_dir()

        # validate input files
        validate_single_file(self.ref_path, '.fa')

        # abort if files exists and --force was not specified
        eprint('Setting up genome reference files in {}'.format(self.out_dir))
        if not delete_or_skip(op.join(self.out_dir, 'CpG.bed.gz'), self.force):
            return

        self.fai_df = self.load_fai()

    def load_fai(self):
        """ Load the fai file to a DataFrame """

        fai_path = self.ref_path + '.fai'

        # If no fai file is found, generate it:
        if not op.isfile(fai_path):
            eprint('fai file not found. Attempting to index {}...'.format(self.ref_path))
            cmd = 'samtools faidx ' + self.ref_path
            subprocess.check_call(cmd, shell=True)
            if op.isfile(fai_path):
                eprint('Generated index file:', fai_path)
            else:
                raise IllegalArgumentError('Failed to generate index file (fai)')

        # Link fa + fai to the output dir
        self.link_file(fai_path, 'genome.fa.fai')
        self.link_file(self.ref_path, 'genome.fa')

        # load fai file
        try:
            df = pd.read_csv(fai_path, sep='\t', header=None, usecols=[0, 1, 2, 4],
                             names=['chr', 'size', 'offset', 'width'])
            # filter invalid chromosomes
            df = df[df.apply(lambda x: is_valid_chrome(x['chr']), axis=1)]

            # sort chromosomes:
            if not self.args.no_sort:
                df = pd.DataFrame(sorted(df['chr'], key=chromosome_order), columns=['chr']).merge(df, how='left')
            return df
        except pd.errors.ParserError as e:
            raise IllegalArgumentError('Invalid fai file.\n{}'.format(e))

    def make_output_dir(self):
        """  equivalent to mkdir -p DIR/references/name/ """

        def mkdir(d):
            if not op.isdir(d):
                os.mkdir(d)
            return d

        return mkdir(op.join(mkdir(op.join(DIR, 'references')), self.name))

    def find_cpgs_loci(self):

        processes = []
        with Pool(self.args.threads) as p:
            for chrom in self.fai_df['chr']:
                params = (chrom, self.ref_path, self.fai_df, self.args.debug)
                processes.append(p.apply_async(load_seq_by_chrom, params))
            p.close()
            p.join()
        df = pd.concat([pr.get() for pr in processes])
        return df

    def run(self):
        eprint('Processing chromosomes...')

        # CpG.bed.gz - Dictionary mapping locus to CpG-Index
        df = self.find_cpgs_loci()
        df.reset_index(drop=True, inplace=True)
        eprint('Composing CpG-Index dictionary..')
        df['site'] = df.index + 1
        self.dump_df(df, 'CpG.bed')
        self.bgzip_tabix_dict(op.join(self.out_dir, 'CpG.bed'))

        # CpG.rev.bin - a reverse dictionary (mapping CpG index to locus
        np.array(df['loc'] + 1, dtype=np.uint32).tofile(op.join(self.out_dir, 'CpG.rev.bin'))

        # chrome.size - number of bases in each chromosome
        self.dump_df(self.fai_df[['chr', 'size']], 'chrome.size')

        # CpG.chrome.size - number of CpGs in each chromosome
        ncgs = self.fai_df[['chr']].copy()
        t = df.groupby('chr')['loc'].nunique().to_frame().reset_index().rename(columns={'loc': 'size'})
        ncgs = ncgs.merge(t, on='chr', how='left').fillna(0)
        ncgs['size'] = ncgs['size'].astype(int)
        self.dump_df(ncgs[['chr', 'size']], 'CpG.chrome.size')

        self.validate_nr_sites(df.shape[0])
        eprint('Finished initialization of genome ', self.name)

    def validate_nr_sites(self, nr_sites):
        if self.args.debug:
            return
        d = {
            'mm9': 13120864,
            'hg19': 28217448
        }
        if self.name in d.keys():
            if nr_sites != d[self.name]:
                eprint('Warning: number of sites of the reference '
                       'genome {} is usually {}, but you got {}'.format(self.name, d[self.name], nr_sites))

    def link_file(self, src, dst):
        if not op.isfile(src):
            raise IllegalArgumentError('Invalid reference genome file: {}'.format(src))
        src = op.abspath(src)
        dst = op.join(self.out_dir, dst)
        cmd = 'ln -s {s} {d}'.format(s=src, d=dst)
        if op.islink(dst):
            os.unlink(dst)
        subprocess.check_output(cmd, shell=True)

    def bgzip_tabix_dict(self, dict_path):
        eprint('bgzip and index...')
        subprocess.check_call('bgzip -@ {} -f '.format(self.args.threads) + dict_path, shell=True)
        subprocess.check_call('tabix -Cf -b 2 -e 2 {}.gz'.format(dict_path), shell=True)

    def dump_df(self, df, path):
        path = op.join(self.out_dir, path)
        df.to_csv(path, index=None, header=False, sep='\t')


def load_seq_by_chrom(chrom, ref_path, fai_df, debug):
    eprint(chrom)

    # get chromosome's location in the fasta
    chrom, size, offset, width = fai_df[fai_df['chr'] == chrom].values[0]

    # load the chromosome's subsequence from fasta
    with open(ref_path, 'r') as f:
        f.seek(offset)
        nr_lines = size // (width - 1) + 1  # number of lines to read for current chromosome
        to_read = nr_lines * width
        if debug:
            to_read = min(to_read, 100 * width)
        txt = f.read(to_read)
    seq = ''.join(s.strip() for s in txt.split('\n')).upper()

    # remove possible trailing characters (belonging to the next chromosome)
    end_pos = seq.rfind('>')
    if end_pos != -1:
        seq = seq[:end_pos]

    # validate sequence length
    if len(seq) != size and not debug:
        raise IllegalArgumentError('Error while loading {} from fasta: '
                                   'read {} bases instead of {}'.format(chrom, len(seq), size))

    # Find CpG sites loci
    tf = pd.DataFrame([m.start() + 1 for m in re.finditer('CG', seq)], columns=['loc'])
    tf['chr'] = chrom
    return tf[['chr', 'loc']]


def chromosome_order(c):
    if not c.startswith('chr'):
        raise IllegalArgumentError('Invalid chromosome' + c)
    c = c[3:]
    if c.isdigit():
        return int(c)
    elif c == 'X':
        return 10000
    elif c == 'Y':
        return 10001
    elif c == 'M':
        return 10002
    return 10003


def is_valid_chrome(chrome):
    # A chromosome is valid if it has the form "chrX", where X is digit(s) or (X,Y,M)
    return bool(re.match(r'^chr([\d]+|[XYM])$', chrome))


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('genome_ref', help='path to a genome *.fa file')
    parser.add_argument('name', help='name of the genome (e.g. hg19, mm9...).\n'
                                     'A directory of this name will be created '
                                     'in references/')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--no_sort', action='store_true',
                        help='If set, keep the chromosome order of the reference genome.\n'
                             'Default behaviour is to sort 1,2,...,10,11,...,X,Y,M')
    parser.add_argument('-@', '--threads', type=int, default=multiprocessing.cpu_count(),
                        help='Number of threads to use (default: multiprocessing.cpu_count)')
    # parser.add_argument('--keep_all', action='store_true', help='Overwrite existing files if existed')
    args = parser.parse_args()
    return args


def main():
    """
    Init genome reference.
    """
    args = parse_args()

    InitGenome(args).run()


if __name__ == '__main__':
    main()
