#!/usr/bin/python3 -u

import os
import os.path as op
import argparse
import numpy as np
import pandas as pd
import subprocess
import shutil
import uuid
import re
from multiprocessing import Pool
from utils_wgbs import IllegalArgumentError, match_maker_tool, patter_tool, \
    add_GR_args, eprint, add_multi_thread_args, EmptyBamError, \
    GenomeRefPaths, validate_single_file, delete_or_skip, check_executable
from init_genome_ref_wgbs import chromosome_order
from pat2beta import pat2beta
from index_wgbs import Indxer
from genomic_region import GenomicRegion

PAT_SUFF = '.pat.gz'

# Minimal Mapping Quality to consider.
# 10 means include only reads w.p. >= 0.9 to be mapped correctly.
# And missing values (255)
MAPQ = 10
FLAGS_FILTER = 1796  # filter flags with these bits
MAX_READ_SIZE = 1000 # extend samtools view region by this size

CHROMS = ['X', 'Y', 'M', 'MT'] + list(range(1, 23))

def extend_region(region, by=MAX_READ_SIZE):
    if ':' not in region:
        return region
    chrom, r = region.split(':')
    start, end = map(int, r.split('-'))
    start = max(1, start - by)
    end += by
    return f'{chrom}:{start}-{end}'

def subprocess_wrap(cmd, debug):
    if debug:
        print(cmd)
        return
    os.system(cmd)
    # p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # output, error = p.communicate()
    # if p.returncode or not output:
        # eprint(cmd)
        # eprint("Failed with subprocess %d\n%s\n%s" % (p.returncode, output.decode(), error.decode()))
        # raise IllegalArgumentError('Failed')

def gen_pat_part(out_path, debug, temp_dir):
    try:
        # if out_path is empty, return None
        if op.getsize(out_path) == 0:
            return
        # sort
        pat_path = out_path + PAT_SUFF
        cmd = f'sort {out_path} -k2,2n -k3,3 '
        if temp_dir:
            cmd += f' -T {temp_dir} '
        cmd += " | uniq -c | awk -v OFS='\t' '{print $2,$3,$4,$1}'"
        cmd += f' | bgzip -f > {pat_path}'
        subprocess_wrap(cmd, debug)

        return pat_path
    except IllegalArgumentError as e:
        return None

def blueprint_legacy(genome, region):
    chrom = region[:region.find(':')]
    # find offset per chrome - how many sites comes before this chromosome
    cf = genome.get_chrom_cpg_size_table()
    cf['size'] = [0] + list(np.cumsum(cf['size']))[:-1]
    chr_offset = cf.set_index('chr').to_dict()['size']

    bppatter_tool = patter_tool.replace('/patter', '/blueprint/patter')
    patter_cmd = f' | {bppatter_tool} {genome.genome_path} {chr_offset[chrom]} '
    patter_cmd += f' --blueprint'
    match_cmd = f' | {match_maker_tool} --drop_singles '
    return patter_cmd, match_cmd


def is_region_empty(view_cmd, verbose):
    # first, if there are no reads in current region, return
    view_cmd += ' | head -1'
    if not subprocess.check_output(view_cmd, shell=True,
            stderr=subprocess.PIPE).decode().strip():
        eprint(f'[wt bam2pat] Skipping region {region}, no reads found')
        if verbose:
            eprint('[wt bam2pat] ' + view_cmd)
        return True
    return False


def proc_chr(bam, out_path, region, genome, chr_offset, paired_end, ex_flags, mapq, debug,
             blueprint, temp_dir, blacklist, whitelist, min_cpg, mbias, verbose):
    """ Convert a temp single chromosome file, extracted from a bam file, into pat """

    # Run patter tool on a single chromosome. out_path will have the following fields:
    # chr   CpG   Pattern   begin_loc   length(bp)


    # use samtools to extract only the reads from 'chrom'
    flag = '-f 3' if paired_end else ''
    view_cmd = f'samtools view {bam} {region} -q {mapq} -F {ex_flags} {flag} '
    if whitelist:
        view_cmd += f' -M -L {whitelist} '
    elif blacklist:
        if not check_executable('bedtools'):
            eprint(f'[wt bam2pat] blacklist flag only works if bedtools is installed')
            raise IllegalArgumentError('Failed')
        view_cmd += f' -b | bedtools intersect -sorted -v -abam stdin -b {blacklist} | samtools view '

    if is_region_empty(view_cmd, verbose):
        return

    # change reads order, s.t paired reads will appear in adjacent lines
    match_cmd = f' | {match_maker_tool} ' if paired_end else ''

    # run patter tool to convert bam reads to pat reads
    patter_cmd = f' | {patter_tool} {genome.dict_path} {extend_region(region)}'
    patter_cmd += f' --min_cpg {min_cpg}'
    if mbias:
        patter_cmd += f' --mbias {out_path}.mb'

    if blueprint:
        patter_cmd, match_cmd = blueprint_legacy(genome, region)
    cmd = view_cmd + match_cmd + patter_cmd + f' > {out_path}'
    if verbose:
        print(cmd)
    subprocess_wrap(cmd, debug)

    return gen_pat_part(out_path, debug, temp_dir)


def validate_bam(bam):

    # validate bam path:
    eprint('[wt bam2pat] bam:', bam)
    if not (op.isfile(bam) and bam.endswith('.bam')):
        eprint(f'[wt bam2pat] Invalid bam: {bam}')
        return False

    # check if bam is sorted by coordinate:
    peek_cmd = f'samtools view -H {bam} | head -1'
    if 'coordinate' not in subprocess.check_output(peek_cmd, shell=True).decode():
        eprint('bam file must be sorted by coordinate')
        return False

    # check if bam is indexed:
    if not (op.isfile(bam + '.bai')):
        eprint('[wt bam2pat] bai file was not found! Generating...')
        if subprocess.call(['samtools', 'index', bam]):
            eprint(f'[wt bam2pat] Failed indexing bam: {bam}')
            return False
    return True


def is_pair_end(bam):
    first_line = subprocess.check_output(f'samtools view {bam} | head -1', shell=True)
    first_line = first_line.decode()
    if len(first_line) == 0:
        raise EmptyBamError('Empty bam file')
    return int(first_line.split('\t')[1]) & 1


class Bam2Pat:
    def __init__(self, args, bam):
        self.args = args
        self.tmp_dir = None
        self.verbose = args.verbose
        self.out_dir = args.out_dir
        self.bam_path = bam
        self.gr = GenomicRegion(args)
        self.start_threads()
        self.cleanup()

    def cleanup(self):
        if self.tmp_dir is not None:
            shutil.rmtree(self.tmp_dir)

    def set_regions(self):
        if self.gr.region_str:
            return [self.gr.region_str]

        cmd = f'samtools idxstats {self.bam_path} | cut -f1 '
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()
        if p.returncode or not output:
            eprint("[wt bam2pat] Failed with samtools idxstats %d\n%s\n%s" % (p.returncode, output.decode(), error.decode()))
            eprint(cmd)
            eprint('[wt bam2pat] falied to find chromosomes')
            return []
        nofilt_chroms = output.decode()[:-1].split('\n')
        filt_chroms = [c for c in nofilt_chroms if 'chr' in c]
        if filt_chroms:
            filt_chroms = [c for c in filt_chroms if re.match(r'^chr([\d]+|[XYM])$', c)]
        else:
            filt_chroms = [c for c in nofilt_chroms if c in CHROMS]
        chroms = list(sorted(filt_chroms, key=chromosome_order))
        if not chroms:
            eprint('[wt bam2pat] Failed retrieving valid chromosome names')
            raise IllegalArgumentError('Failed')

        return chroms

    def set_lists(self):
        # black/white lists:
        blacklist = self.args.blacklist
        whitelist = self.args.whitelist
        if blacklist == True:
            blacklist = GenomeRefPaths(self.args.genome).blacklist
        elif whitelist == True:
            whitelist = GenomeRefPaths(self.args.genome).whitelist
        if blacklist:
            validate_single_file(blacklist)
        elif whitelist:
            validate_single_file(whitelist)
        if self.verbose:
            eprint(f'[wt bam2pat] blacklist: {blacklist}')
            eprint(f'[wt bam2pat] whitelist: {whitelist}')
        return blacklist, whitelist

    def start_threads(self):
        """ Parse each chromosome file in a different process,
            and concatenate outputs to pat files """

        blist, wlist = self.set_lists()
        name = op.join(self.out_dir, op.basename(self.bam_path)[:-4])
        # build temp dir:
        name = op.splitext(op.basename(self.bam_path))[0]
        self.tmp_dir = op.join(self.out_dir,
                f'{name}.{str(uuid.uuid4())[:8]}.PID{os.getpid()}')
        os.mkdir(self.tmp_dir)
        tmp_prefix = op.join(self.tmp_dir, name)

        # find offset per chrome - how many sites comes before this chromosome
        cf = GenomeRefPaths(self.gr.genome_name).get_chrom_cpg_size_table()
        cf['size'] = [0] + list(np.cumsum(cf['size']))[:-1]
        chr_offset = cf.set_index('chr').to_dict()['size']

        params = []
        cur_regions = self.set_regions()
        try:
            for c in cur_regions:
                out_path = f'{tmp_prefix}.{c}.out'
                par = (self.bam_path, out_path, c, self.gr.genome, chr_offset,
                       is_pair_end(self.bam_path), self.args.exclude_flags,
                       self.args.mapq, self.args.debug, self.args.blueprint,
                       self.args.temp_dir, blist, wlist, self.args.min_cpg,
                       self.args.mbias, self.verbose)
                params.append(par)

            if len(cur_regions) == 1 and self.args.threads == 1:
                cur_params = params[0]
                pat_parts = [proc_chr(*cur_params)]
            else:
                if not params:
                    raise IllegalArgumentError('Empty bam file')

                p = Pool(self.args.threads)
                pat_parts = p.starmap(proc_chr, params)
                p.close()
                p.join()
        except EmptyBamError as e:
            if len(cur_regions) == 1 and self.args.threads == 1:
                eprint("")
                pat_parts = []
            else:
                raise e

        self.mbias_merge(name, pat_parts)
        self.concat_parts(name, pat_parts)

    def validate_parts(self, pat_parts):
        # validate parts:
        pat_parts = [p for p in pat_parts if p]  # in case only some of the parts are empty

        for part in pat_parts:
            if not part:            # subprocess threw exception or no reads found
                eprint('[wt bam2pat] threads failed')
                return False
            if op.getsize(part) == 0:   # empty pat was created
                eprint('[wt bam2pat] threads failed')
                return False
        if not ''.join(pat_parts):      # all parts are empty
            eprint('[wt bam2pat] No reads found in bam file. No pat file is generated')
            return False
        return True

    def mbias_merge(self, name, pat_parts):
        if not self.args.mbias:
            return
        try:
            mdir = op.join(self.out_dir, name) + '.mbias'
            if not op.isdir(mdir):
                os.mkdir(mdir)
            tpaths = []
            for x in ['OB', 'OT']:
                mbias_parts = [p.replace('.pat.gz', f'.mb.{x}.txt') for p in pat_parts if p]
                mbias_parts = [pd.read_csv(m, sep='\t') for m in mbias_parts]
                df = mbias_parts[0]
                for m in mbias_parts[1:]:
                    df += m
                cpath = op.join(mdir, name) + f'.mbias.{x}.txt'
                df.to_csv(cpath, sep='\t', index=None)
                tpaths.append(cpath)
            from mbias_plot import plot_mbias
            plot_mbias(tpaths, mdir)
        except Exception as e:
            eprint('[wt bam2pat] failed in mbias')
            eprint(e)


    def concat_parts(self, name, pat_parts):

        if not self.validate_parts(pat_parts):
            return

        # Concatenate chromosome files
        pat_path = op.join(self.out_dir, name) + PAT_SUFF
        os.system('cat ' + ' '.join(pat_parts) + ' > ' + pat_path)

        if not op.isfile(pat_path):
            eprint(f'[wt bam2pat] failed to generate {pat_path}')
            return

        # index the pat file (.csi)
        eprint('[wt bam2pat] indexing and generating beta file...')
        Indxer(pat_path, threads=self.args.threads).run()
        eprint(f'[wt bam2pat] generated {pat_path}')

        # generate beta file
        if not self.args.no_beta:
            beta_path = pat2beta(f'{pat_path}', self.out_dir, args=self.args)
            eprint(f'[wt bam2pat] generated {beta_path}')


def parse_bam2pat_args(parser):
    parser.add_argument('--no_beta', action='store_true', help='Do not generate a beta fils, only pat')
    parser.add_argument('-l', '--lbeta', action='store_true', help='Use lbeta file (uint16) instead of beta (uint8)')
    parser.add_argument('-T', '--temp_dir', help='passed to unix sort. Useful in case bam file is very large')
    lists = parser.add_mutually_exclusive_group()
    lists.add_argument('--blacklist', help='bed file. Ignore reads overlapping this bed file',
                        nargs='?', const=True, default=False)
    lists.add_argument('-L', '--whitelist', help='bed file. Consider only reads overlapping this bed file',
                        nargs='?', const=True, default=False)
    parser.add_argument('--mbias', '-mb', action='store_true',
            help='Output mbias plots. Only paired-end data is supported')
    parser.add_argument('--blueprint', '-bp', action='store_true',  # TODO put it back
            help='filter bad bisulfite conversion reads if <90 percent of CHs are converted')


def add_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('bam', nargs='+')
    add_GR_args(parser)
    parser.add_argument('--out_dir', '-o', default='.')
    parser.add_argument('--min_cpg', type=int, default=1,
                help='Reads covering less than MIN_CPG sites are removed [1]')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--force', '-f', action='store_true', help='overwrite existing files if exists')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('-F', '--exclude_flags', type=int,
                        help='flags to exclude from bam file (samtools view parameter) ' \
                             f'[{FLAGS_FILTER}]', default=FLAGS_FILTER)
    parser.add_argument('-q', '--mapq', type=int,
                        help=f'Minimal mapping quality (samtools view parameter) [{MAPQ}]',
                        default=MAPQ)
    add_multi_thread_args(parser)

    return parser


def parse_args(parser):
    parse_bam2pat_args(parser)
    args = parser.parse_args()
    return args


def main():
    """
    Run the WGBS pipeline to generate pat & beta files out of an input bam file
    """
    parser = add_args()
    args = parse_args(parser)
    # validate output dir:
    if not op.isdir(args.out_dir):
        raise IllegalArgumentError(f'Invalid output dir: {args.out_dir}')

    for bam in args.bam:
        if not validate_bam(bam):
            eprint(f'[wt bam2pat] Skipping {bam}')
            continue

        pat = op.join(args.out_dir, op.basename(bam)[:-4] + PAT_SUFF)
        if not delete_or_skip(pat, args.force):
            continue
        Bam2Pat(args, bam)


if __name__ == '__main__':
    main()
