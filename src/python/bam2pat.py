#!/usr/bin/python3 -u

import os
import os.path as op
import subprocess
import shutil
import uuid
import argparse
from multiprocessing import Pool
import numpy as np
import pandas as pd
from utils_wgbs import IllegalArgumentError, match_maker_tool, patter_tool, \
    add_GR_args, eprint, add_multi_thread_args, EmptyBamError, \
    validate_single_file, delete_or_skip, check_executable, \
    add_no_beta_arg, validate_local_exe, mkdirp, pretty_name
from init_genome import chromosome_order
from pat2beta import pat2beta
from index import Indxer
from genomic_region import GenomicRegion

PAT_SUFF = '.pat.gz'

# Minimal Mapping Quality to consider.
# 10 means include only reads w.p. >= 0.9 to be mapped correctly.
# And missing values (255)
MAPQ = 10
FLAGS_FILTER = 1796  # filter flags with these bits
MAX_READ_SIZE = 1000 # extend samtools view region by this size

# currently we are missing the far mates
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
        eprint(cmd)
        return
    os.system(cmd)
    # p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # output, error = p.communicate()
    # if p.returncode or not output:
        # eprint(cmd)
        # eprint("Failed with subprocess %d\n%s\n%s" % (p.returncode, output.decode(), error.decode()))
        # raise IllegalArgumentError('Failed')


def set_regions(bam_path, gr, tmp_dir=None):
    # if user specified a region, just use it
    if gr.region_str:
        return [gr.region_str]

    # get all chromosomes from the reference genome:
    ref_chroms = gr.genome.get_chroms()
    if bam_path.endswith('.cram'):
        return list(sorted(ref_chroms, key=chromosome_order))

    # get all chromosomes present in the bam file header
    cmd = f'samtools idxstats {bam_path} | cut -f1 '
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.returncode or not output:
        eprint("[wt bam2pat] Failed with samtools idxstats %d\n%s\n%s" % (p.returncode, output.decode(), error.decode()))
        eprint(cmd)
        eprint('[wt bam2pat] falied to find chromosomes')
        return []
    bam_chroms = output.decode()[:-1].split('\n')

    # intersect the chromosomes from the bam and from the reference
    intersected_chroms = list(set(bam_chroms) & set(ref_chroms))
    if not intersected_chroms:
        msg = '[wt bam2pat] Failed retrieving valid chromosome names. '
        msg += 'Perhaps you are using a wrong genome reference. '
        msg += 'Try running:\n\t\twgbstools set_default_ref -ls'
        eprint(msg)
        tmpdir_cleanup(tmp_dir)
        raise IllegalArgumentError('Failed')

    return list(sorted(intersected_chroms, key=chromosome_order))


def tmpdir_cleanup(tmp_dir):
    if tmp_dir is not None:
        shutil.rmtree(tmp_dir)


def gen_pat_part(out_path, debug, temp_dir):
    try:
        # if out_path is empty or missing, return None
        if not op.isfile(out_path):
            eprint(f'[wt bam2pat] Failed in patter: file {out_path} is missing')
            return
        if op.getsize(out_path) == 0:
            eprint(f'[wt bam2pat] Failed in patter: file {out_path} is empty')
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
    except IllegalArgumentError:
        return None

def blueprint_legacy(genome, region, paired_end):
    if not op.isfile(genome.genome_path):
        eprint(f'[ wt bam2pat ] Error: not genome reference fasta file: {genome.genome_path}')
        raise IllegalArgumentError('Failed')
    chrom = region
    if ':' in region:
        chrom = region[:region.find(':')]
    # find offset per chrome - how many sites comes before this chromosome
    cf = genome.get_chrom_cpg_size_table()
    cf['size'] = [0] + list(np.cumsum(cf['size']))[:-1]
    chr_offset = cf.set_index('chr').to_dict()['size']

    bppatter_tool = patter_tool.replace('/patter', '/blueprint/patter')
    patter_cmd = f' | {bppatter_tool} {genome.genome_path} {chr_offset[chrom]} '
    patter_cmd += ' --blueprint'
    match_cmd = f' | {match_maker_tool} --drop_singles ' if paired_end else ''
    return patter_cmd, match_cmd


def is_region_empty(view_cmd, region, verbose):
    # check if there are reads in the bam file for the requested region
    view_cmd += ' | head -1'
    if not subprocess.check_output(view_cmd, shell=True,
            stderr=subprocess.PIPE).decode().strip():
        if verbose:
            eprint(f'[wt bam2pat] Skipping region {region}, no reads found')
            eprint('[wt bam2pat] ' + view_cmd)
        return True
    return False


def proc_chr(bam, out_path, region, genome, paired_end, ex_flags, in_flags, mapq, debug,
             blueprint, clip, temp_dir, blacklist, whitelist, min_cpg, mbias, nanopore,
             np_thresh, verbose):
    """ Convert a temp single chromosome file, extracted from a bam file, into pat """

    # Run patter tool on a single chromosome (or region). out_path will have the following fields:
    # chr   CpG   Pattern   begin_loc   length(bp)


    # use samtools to extract only the reads from 'chrom'
    if in_flags is None:
        in_flags = '-f 3' if paired_end else ''
    else:
        in_flags = f'-f {in_flags}'

    view_cmd = f'samtools view {bam} {region} -q {mapq} -F {ex_flags} {in_flags} -T {genome.genome_path}'
    if whitelist:
        view_cmd += f' -M -L {whitelist} '
    elif blacklist:
        if not check_executable('bedtools'):
            eprint('[wt bam2pat] blacklist flag only works if bedtools is installed')
            raise IllegalArgumentError('Failed')
        view_cmd += f' -b | bedtools intersect -sorted -v -abam stdin -b {blacklist} | samtools view '

    # first, if there are no reads in current region, return
    if is_region_empty(view_cmd, region, verbose):
        return

    # change reads order, s.t paired reads will appear in adjacent lines
    match_cmd = f' | {match_maker_tool} ' if paired_end else ''

    # run patter tool to convert bam reads to pat reads
    patter_cmd = f' | {patter_tool} {genome.dict_path} {extend_region(region)}'
    patter_cmd += f' --min_cpg {min_cpg} --clip {clip}'
    if mbias:
        patter_cmd += f' --mbias {out_path}.mb'
    if nanopore:
        patter_cmd += f' --nanopore --np_thresh {np_thresh} '

    if blueprint:
        patter_cmd, match_cmd = blueprint_legacy(genome, region, paired_end)
    cmd = view_cmd + match_cmd + patter_cmd + f' > {out_path}'
    if verbose:
        print(cmd)
    subprocess_wrap(cmd, debug)

    return gen_pat_part(out_path, debug, temp_dir)


def validate_bam(bam):

    # validate bam path:
    eprint('[wt bam2pat] bam:', bam)
    if not (op.isfile(bam) and bam.endswith(('.bam', '.cram'))):
        eprint(f'[wt bam2pat] Invalid bam: {bam}')
        return False
    return True


def is_bam_sorted(bam):

    # check if bam is sorted by coordinate:
    peek_cmd = f'samtools view -H {bam}| head -1'
    hd_line = subprocess.check_output(peek_cmd, shell=True).decode()
    if hd_line.startswith('@HD') and 'coordinate' not in hd_line:
        eprint(f'[wt bam2pat] WARNING: based on the @HD, bam file is not sorted: {bam}')
        return False

    # check if bam is indexed:
    is_indexed = (op.isfile(bam + '.bai') or op.isfile(bam + '.csi') or op.isfile(bam + '.crai'))
    if not is_indexed:
        eprint('[wt bam2pat] WARNING: index file (bai/csi) not found! Attempting to generate bai...')
        if subprocess.call(['samtools', 'index', bam]):
            eprint(f'[wt bam2pat] Failed indexing bam: {bam}')
            eprint('              Make sure the bam file is sorted and indexed')
            return False

    return True


def is_pair_end(bam, genome):
    first_line = subprocess.check_output(f'samtools view {bam} -T {genome.genome_path} | head -1', shell=True)
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
        self.PE = None
        self.start_threads()
        tmpdir_cleanup(self.tmp_dir)

    def set_lists(self):
        # black/white lists:
        blacklist = self.args.blacklist
        whitelist = self.args.whitelist
        if blacklist is True:
            blacklist = self.gr.genome.blacklist
        elif whitelist is True:
            whitelist = self.gr.genome.whitelist
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

        self.PE = is_pair_end(self.bam_path, self.gr.genome)
        blist, wlist = self.set_lists()
        # build temp dir:
        name = pretty_name(self.bam_path)
        self.tmp_dir = op.join(self.out_dir,
                f'{name}.PID{os.getpid()}.{str(uuid.uuid4())[:8]}')
        os.mkdir(self.tmp_dir)
        tmp_prefix = op.join(self.tmp_dir, name)

        params = []
        cur_regions = set_regions(self.bam_path, self.gr, self.tmp_dir)
        try:
            for c in cur_regions:
                out_path = f'{tmp_prefix}.{c}.out'
                par = (self.bam_path, out_path, c, self.gr.genome,
                       self.PE, self.args.exclude_flags,
                       self.args.include_flags,
                       self.args.mapq, self.args.debug, self.args.blueprint, self.args.clip,
                       self.args.temp_dir, blist, wlist, self.args.min_cpg,
                       self.args.mbias, self.args.nanopore, self.args.np_thresh, self.verbose)
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
                tmpdir_cleanup(self.tmp_dir)
                raise e

        self.mbias_merge(name, pat_parts)
        self.concat_parts(name, pat_parts)

    @staticmethod
    def validate_parts(pat_parts):
        # validate parts:
        pat_parts = [p for p in pat_parts if p]  # in case only some of the parts are empty

        for part in pat_parts:
            if not part:            # subprocess threw exception or no reads found
                eprint('[wt bam2pat] threads failed')
                return []
            if op.getsize(part) == 0:   # empty pat was created
                eprint('[wt bam2pat] threads failed')
                return []
        if not ''.join(pat_parts):      # all parts are empty
            eprint('[wt bam2pat] No reads found. No pat file is generated')
            return []
        return pat_parts

    def mbias_merge(self, name, pat_parts):
        if not self.args.mbias:
            return
        try:
            mdir = op.join(self.out_dir, name) + '.mbias'
            mkdirp(mdir)
            tpaths = []
            for x in ['OT', 'OB']:
                mbias_parts = [p.replace('.pat.gz', f'.mb.{x}.txt') for p in pat_parts if p]
                mbias_parts = [pd.read_csv(m, sep='\t') for m in mbias_parts]
                df = mbias_parts[0]
                for m in mbias_parts[1:]:
                    df += m
                cpath = op.join(mdir, name) + f'.mbias.{x}.txt'
                df.to_csv(cpath, sep='\t', index=None)
                tpaths.append(cpath)
            from mbias_plot import plot_mbias
            plot_mbias(tpaths, mdir, self.PE)
        except Exception as e:
            eprint('[wt bam2pat] failed in mbias')
            eprint(e)


    def concat_parts(self, name, pat_parts):

        pat_parts = self.validate_parts(pat_parts)
        if not pat_parts:
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
    add_no_beta_arg(parser)
    parser.add_argument('-l', '--lbeta', action='store_true', help='Use lbeta file (uint16) instead of beta (uint8)')
    parser.add_argument('-T', '--temp_dir', help='passed to unix sort. Useful in case bam file is very large')
    lists = parser.add_mutually_exclusive_group()
    lists.add_argument('--blacklist', help='bed file. Ignore reads overlapping this bed file',
                        nargs='?', const=True, default=False)
    lists.add_argument('-L', '--whitelist', help='bed file. Consider only reads overlapping this bed file',
                        nargs='?', const=True, default=False)
    parser.add_argument('--mbias', '-mb', action='store_true',
            help='Output mbias plots. Only paired-end data is supported')
    parser.add_argument('--blueprint', '-bp', action='store_true',
            help='filter bad bisulfite conversion reads if <90 percent of CHs are converted')
    parser.add_argument('--nanopore', '-np', action='store_true',
            help='BETA VERSION: Input bam is of Oxford Nanopore format. Call methylation from MM & ML fields.')
    parser.add_argument('--np_thresh', type=float, default=0.67,
                        help='For Nanopore format: probability cutoff, between 0 to 1. [0.67]')


def add_samtools_view_flags(parser):
    parser.add_argument('--include_flags', type=int,
                        help='flags to include from bam file (samtools view parameter -f) ' \
                             '[3 for PE, None for SE]')
    parser.add_argument('-F', '--exclude_flags', type=int,
                        help='flags to exclude from bam file (samtools view parameter -F) ' \
                             f'[{FLAGS_FILTER}]', default=FLAGS_FILTER)
    parser.add_argument('-q', '--mapq', type=int,
                        help=f'Minimal mapping quality (samtools view parameter) [{MAPQ}]',
                        default=MAPQ)


def add_args(parser):
    parser.add_argument('bam', nargs='+')
    add_GR_args(parser)
    add_samtools_view_flags(parser)
    parser.add_argument('--out_dir', '-o', default='.')
    parser.add_argument('--min_cpg', type=int, default=1,
                help='Reads covering less than MIN_CPG sites are removed [1]')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--force', '-f', action='store_true', help='overwrite existing files if exists')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--clip', type=int, default=0,
                        help='Clip for each read the first and last CLIP characters [0]')
    add_multi_thread_args(parser)

    return parser


def parse_args(parser):
    parse_bam2pat_args(parser)
    args = parser.parse_args()
    return args


def validate_np_thresh(args):
    if 'np_thresh' in args:
        if not 0 < args.np_thresh < 1:
            raise IllegalArgumentError('Invalid np_thresh range: must be in range (0,1)')


def main():
    """
    Run the WGBS pipeline to generate pat & beta files out of an input bam file
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser = add_args(parser)
    args = parse_args(parser)
    # validate output dir:
    if not op.isdir(args.out_dir):
        raise IllegalArgumentError(f'Invalid output dir: {args.out_dir}')

    validate_local_exe(match_maker_tool)
    validate_local_exe(patter_tool)
    validate_np_thresh(args)

    for bam in args.bam:
        if not (validate_bam(bam) and is_bam_sorted(bam)):
            eprint(f'[wt bam2pat] Skipping {bam}')
            continue

        pat = op.join(args.out_dir, pretty_name(bam) + PAT_SUFF)
        if not delete_or_skip(pat, args.force):
            continue
        Bam2Pat(args, bam)


if __name__ == '__main__':
    main()
