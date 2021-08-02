#!/usr/bin/python3 -u

import argparse
import re
import shutil
import numpy as np
import pandas as pd
import os.path as op
import os
from itertools import groupby
import subprocess
from multiprocessing import Pool
import multiprocessing
from utils_wgbs import validate_single_file, eprint, IllegalArgumentError, DIR, add_multi_thread_args


class InitGenome:
    def __init__(self, args):
        self.args = args
        self.ref_path = args.genome_ref
        self.force = args.force
        self.name = args.name

        # validate input files
        validate_single_file(self.ref_path, '.fa')  # todo: support bgzipped reference FASTA
                                                    # It's half implemented - only patter does not support it
        self.out_dir = self.setup_dir()
        self.fai_df = self.load_fai()

    def setup_dir(self):
        out_dir = op.join(op.join(op.join(DIR, 'references')), self.name)
        # abort if files exists and --force was not specified
        if op.isdir(out_dir):
            if not self.force:
                msg = f'[wt init] Error: genome {self.name} already exists ({out_dir}). Use -f to overwrite it.'
                raise IllegalArgumentError(msg)
            else:
                shutil.rmtree(out_dir)
        os.makedirs(out_dir, exist_ok=True)
        eprint(f'[wt init] Setting up genome reference files in {out_dir}')
        return out_dir

    def load_fai(self):
        """ Load the fai file to a DataFrame """

        fai_path = self.ref_path + '.fai'

        # If no fai file is found, generate it:
        if not op.isfile(fai_path):
            eprint(f'[wt init] fai file not found. Attempting to index {self.ref_path}')
            cmd = f'samtools faidx {self.ref_path}'
            output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode()
            if self.ref_path.endswith('.gz') and 'please use bgzip' in output:
                msg = f'[wt init] Seems like your reference FASTA cannot be indexed with samtools faidx.\n' \
                        ' Try one of the following:\n' \
                        ' 1. decompress it (gunzip {self.ref_path}) and try again' \
                        ' 2. change the compression to bgzip:\n' \
                        'gzip {self.ref_path} && bgzip {self.ref_path[:-3]}'
                eprint(msg)
                raise IllegalArgumentError('[wt init] Invalid reference FASTA')
            if op.isfile(fai_path):
                eprint(f'[wt init] Generated index file: {fai_path}')
            else:
                raise IllegalArgumentError('[wt init] Failed to generate index file (fai)')

        # Link fa + fai to the output dir
        # fasta_name = 'genome.fa' + '.gz' if self.ref_path.endswith('.gz') else ''
        fasta_name = 'genome.fa'
        self.link_file(self.ref_path, fasta_name)
        self.link_file(fai_path, fasta_name + '.fai')

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
            raise IllegalArgumentError(f'Invalid fai file.\n{e}')

    def find_cpgs_loci(self):
        p = Pool(self.args.threads)
        params = [(self.ref_path, c) for c in self.fai_df['chr']]
        arr = p.starmap(load_seq_by_chrom, params)
        p.close()
        p.join()
        return pd.concat(arr)

    def run(self):
        eprint('[wt init] Processing chromosomes...')

        # CpG.bed.gz - Dictionary mapping locus to CpG-Index
        df = self.find_cpgs_loci()
        df.reset_index(drop=True, inplace=True)
        eprint('[wt init] Building CpG-Index dictionary...')
        df['site'] = df.index + 1
        self.dump_df(df, 'CpG.bed')
        cpg_dict = self.bgzip_tabix_dict(op.join(self.out_dir, 'CpG.bed'))

        # rev.CpG.bed.gz - symbolic link to the dictionary, 
        # with an index file corresponds to the 3rd column, the CpG-Index,
        # To map CpG to locus
        rev_dict = op.join(self.out_dir, 'rev.CpG.bed.gz')
        os.symlink(cpg_dict, rev_dict)
        subprocess.check_call(f'tabix -b 3 -e 3 {rev_dict}', shell=True)

        # chrome.size - number of bases in each chromosome
        self.dump_df(self.fai_df[['chr', 'size']], 'chrome.size')

        # CpG.chrome.size - number of CpGs in each chromosome
        ncgs = self.fai_df[['chr']].copy()
        t = df.groupby('chr')['loc'].nunique().to_frame().reset_index().rename(columns={'loc': 'size'})
        ncgs = ncgs.merge(t, on='chr', how='left').fillna(0)
        ncgs['size'] = ncgs['size'].astype(int)
        self.dump_df(ncgs[['chr', 'size']], 'CpG.chrome.size')

        self.validate_nr_sites(df.shape[0])
        eprint(f'[wt init] Finished initialization of genome {self.name}')

    def validate_nr_sites(self, nr_sites):
        if self.args.debug:
            return
        d = {
            'mm9': 13120864,
            'hg19': 28217448
        }
        if self.name in d.keys():
            if nr_sites != d[self.name]:
                msg = f'[wt init] WARNING: number of sites of the reference genome '
                msg += f'{self.name} is usually {d[self.name]}, but you got {nr_sites}'
                eprint(msg)

    def link_file(self, src, dst):
        if not op.isfile(src):
            raise IllegalArgumentError(f'[wt init] Invalid reference genome file: {src}')
        src = op.abspath(src)
        dst = op.join(self.out_dir, dst)
        cmd = f'ln -s {src} {dst}'
        if op.islink(dst):
            os.unlink(dst)
        subprocess.check_call(cmd, shell=True)

    def bgzip_tabix_dict(self, dict_path):
        eprint('[wt init] bgzip and index...')
        subprocess.check_call(f'bgzip -@ {self.args.threads} -f {dict_path}', shell=True)
        subprocess.check_call(f'tabix -Cf -b 2 -e 2 {dict_path}.gz', shell=True)
        return dict_path + '.gz'

    def dump_df(self, df, path):
        path = op.join(self.out_dir, path)
        df.to_csv(path, index=None, header=False, sep='\t')


def load_seq_by_chrom(ref_path, chrom):
    eprint(f'[wt init] {chrom}')

    # load the chromosome's subsequence from fasta
    cmd = f'samtools faidx {ref_path} {chrom} | tail -n +2'
    txt = subprocess.check_output(cmd, shell=True).decode()
    seq = ''.join(s.strip() for s in txt.split('\n')).upper()

    # Find CpG sites loci
    tf = pd.DataFrame([m.start() + 1 for m in re.finditer('CG', seq)], columns=['loc'])
    tf['chr'] = chrom
    return tf[['chr', 'loc']]


def chromosome_order(c):
    if c.startswith('chr'):
        # raise IllegalArgumentError('Invalid chromosome' + c)
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
    parser.add_argument('genome_ref', help='path to a genome FASTA file')
    parser.add_argument('name', help='name of the genome (e.g. hg19, mm9...).\n'
                                     'A directory of this name will be created '
                                     'in references/')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing files if existed')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--no_sort', action='store_true',
                        help='If set, keep the chromosome order of the reference genome.\n'
                             'Default behaviour is to sort 1,2,...,10,11,...,X,Y,M')
    add_multi_thread_args(parser)
    # parser.add_argument('--keep_all', action='store_true', help='Overwrite existing files if existed')
    args = parser.parse_args()
    return args


def main():
    """
    Init genome reference.
    Note: we currently support only chromosomes starting with "chr". I.e. "chr1" and not "1".
    """
    args = parse_args()

    InitGenome(args).run()


if __name__ == '__main__':
    main()
