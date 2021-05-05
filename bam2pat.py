#!/usr/bin/python3 -u


import os
import os.path as op
import argparse
import subprocess
import shutil
import uuid
import re
from multiprocessing import Pool
from utils_wgbs import IllegalArgumentError, match_maker_tool, patter_tool, add_GR_args, eprint, \
        add_multi_thread_args, mult_safe_remove, collapse_pat_script, GenomeRefPaths, validate_single_file
from init_genome_ref_wgbs import chromosome_order
from pat2beta import pat2beta
from genomic_region import GenomicRegion

PAT_SUFF = '.pat'
UNQ_SUFF = '.unq'

# Minimal Mapping Quality to consider.
# 10 means include only reads w.p. >= 0.9 to be mapped correctly.
# And missing values (255)
MAPQ = 10
FLAGS_FILTER = 1796  # filter flags with these bits

# todo: unsorted / sorted by name
CHROMS = ['X', 'Y', 'M', 'MT'] + list(range(1, 23))


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

def pat(out_path, debug, temp_dir):
    try:
        # sort
        pat_path = out_path + PAT_SUFF
        cmd = f'sort {out_path} -k2,2n -k3,3 '
        if temp_dir:
            cmd += f' -T {temp_dir} '
        cmd += f' | sed -e "s/$/\t1/" | {collapse_pat_script} - > {pat_path}'
        subprocess_wrap(cmd, debug)

        mult_safe_remove([out_path])

        return pat_path, None
    except IllegalArgumentError as e:
        return None

def pat_unq(out_path, debug, unq, temp_dir):
    try:
        # sort
        cmd = f'sort {out_path} -k2,2n -k3,3 -o {out_path}'
        if temp_dir:
            cmd += f' -T {temp_dir} '
        subprocess_wrap(cmd, debug)

        # break output file into pat and unq:
        # pat file:
        pat_path = out_path + PAT_SUFF
        cmd = 'awk \'BEGIN{OFS="\t"}{print $1,$2,$3,1}\' ' + f'{out_path} | {collapse_pat_script} - > {pat_path}'
        subprocess_wrap(cmd, debug)

        # unq file:
        unq_path = out_path + UNQ_SUFF
        if unq:
            cmd = "sort {} -k4,4n -k3,3 -o {}".format(out_path, out_path)
            if temp_dir:
                cmd += ' -T {} '.format(temp_dir)
            subprocess_wrap(cmd, debug)
            cmd = 'awk \'{print $1,$4,$5,$3}\' ' + out_path + ' | uniq -c | awk \'{OFS="\\t"; print $2,$3,$4,$5,$1}\' > ' + \
                  unq_path
            subprocess_wrap(cmd, debug)

        mult_safe_remove([out_path])

        return pat_path, unq_path
    except IllegalArgumentError as e:
        return None


def proc_chr(bam, out_path, region, genome, paired_end, ex_flags, mapq, debug,
             unq, blueprint, temp_dir, blacklist, whitelist, verbose):
    """ Convert a temp single chromosome file, extracted from a bam file,
        into two output files: pat and unq."""

    # Run patter tool on a single chromosome. out_path will have the following fields:
    # chr   CpG   Pattern   begin_loc   length(bp)


    # use samtools to extract only the reads from 'chrom'
    flag = '-f 3' if paired_end else ''
    cmd = f'samtools view {bam} {region} -q {mapq} -F {ex_flags} {flag} '
    if whitelist:
        cmd += f' -M -L {whitelist} '
    elif blacklist:
        cmd += f' -b | bedtools intersect -sorted -v -abam stdin -b {blacklist} | samtools view '
    if debug:
        cmd += ' | head -200 '
    if paired_end:
        # change reads order, s.t paired reads will appear in adjacent lines
        cmd += f' | {match_maker_tool} '

    # first, if there are no reads in current region, return
    validation_cmd = cmd + ' | head -1'
    if not subprocess.check_output(validation_cmd, shell=True, stderr=subprocess.PIPE).decode().strip():
        eprint(f'[bam2pat] Skipping region {region}, no reads found')
        if verbose:
            eprint('[bam2pat] ' + validation_cmd)
        return '', ''

    cmd += f' | {patter_tool} {genome.genome_path} {genome.chrom_cpg_sizes} '
    if unq:
        cmd += ' --unq '
    if blueprint:
        cmd += ' --blueprint '
    cmd += f' > {out_path}'
    if verbose:
        print(cmd)
    subprocess_wrap(cmd, debug)

    if unq:
        return pat_unq(out_path, debug, unq, temp_dir)
    else:
        return pat(out_path, debug, temp_dir)


class Bam2Pat:
    def __init__(self, args):
        self.args = args
        self.tmp_dir = None
        self.verbose = args.verbose
        self.out_dir = args.out_dir
        self.bam_path = args.bam_path
        self.debug = args.debug
        self.gr = GenomicRegion(args)
        self.validate_input()
        self.start_threads()
        self.cleanup()

    def cleanup(self):
        if self.tmp_dir is not None:
            shutil.rmtree(self.tmp_dir)

    def validate_input(self):

        # validate bam path:
        eprint('[bam2pat] bam:', self.bam_path)
        if not (op.isfile(self.bam_path) and self.bam_path.endswith('.bam')):
            raise IllegalArgumentError('[bam2pat] Invalid bam: {}'.format(self.bam_path))

        # check if bam is sorted by coordinate:
        peek_cmd = f'samtools view -H {self.bam_path} | head -1'
        so = subprocess.PIPE
        if 'coordinate' not in subprocess.check_output(peek_cmd, shell=True).decode():
            raise IllegalArgumentError('bam file must be sorted by coordinate')

        # check if bam is indexed:
        if not (op.isfile(self.bam_path + '.bai')):
            eprint('[bam2pat] bai file was not found! Generating...')
            r = subprocess.call(['samtools', 'index', self.bam_path])
            if r:
                raise IllegalArgumentError(f'Failed indexing bam: {self.bam_path}')

        # validate output dir:
        if not (op.isdir(self.out_dir)):
            raise IllegalArgumentError(f'Invalid output dir: {self.out_dir}')

    def is_pair_end(self):
        first_line = subprocess.check_output(f'samtools view {self.bam_path} | head -1', shell=True)
        return int(first_line.decode().split('\t')[1]) & 1

    def set_regions(self):
        if self.gr.region_str:
            return [self.gr.region_str]
        else:
            cmd = f'samtools idxstats {self.bam_path} | cut -f1 '
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = p.communicate()
            if p.returncode or not output:
                eprint(cmd)
                eprint("[bam2pat] Failed with samtools idxstats %d\n%s\n%s" % (p.returncode, output.decode(), error.decode()))
                eprint('[bam2pat] falied to find chromosomes')
                return []
            nofilt_chroms = output.decode()[:-1].split('\n')
            filt_chroms = [c for c in nofilt_chroms if 'chr' in c]
            if not filt_chroms:
                filt_chroms = [c for c in nofilt_chroms if c in CHROMS]
            else:
                filt_chroms = [c for c in filt_chroms if re.match(r'^chr([\d]+|[XYM])$', c)]
            chroms = list(sorted(filt_chroms, key=chromosome_order))
            if not chroms:
                eprint('[bam2pat] Failed retrieving valid chromosome names')
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
            and concatenate outputs to pat and unq files """

        blist, wlist = self.set_lists()
        name = op.join(self.out_dir, op.basename(self.bam_path)[:-4])
        # build temp dir:
        name = op.splitext(op.basename(self.bam_path))[0]
        uname = f'{name}.{str(uuid.uuid4())[:8]}.PID{os.getpid()}'
        self.tmp_dir = op.join(self.args.out_dir, uname)
        os.mkdir(self.tmp_dir)
        tmp_prefix = op.join(self.tmp_dir, name)

        processes = []
        nr_threads = self.args.threads  # todo smarter default!
        with Pool(nr_threads) as p:
            for c in self.set_regions():
                out_path = f'{tmp_prefix}.{c}.out'
                params = (self.bam_path, out_path, c, self.gr.genome,
                          self.is_pair_end(), self.args.exclude_flags,
                          self.args.mapq, self.debug, self.args.unq, self.args.blueprint,
                          self.args.temp_dir, blist, wlist, self.verbose)
                processes.append(p.apply_async(proc_chr, params))
            if not processes:
                raise IllegalArgumentError('Empty bam file')
            p.close()
            p.join()
        res = [pr.get() for pr in processes]  # [(pat_path, unq_path) for each chromosome]
        if None in res:
            eprint('[bam2pat] threads failed')
            return
        if not ''.join(p for p, u in res):
            eprint('[bam2pat] No reads found in bam file. No pat file is generated')
            return

        # Concatenate chromosome files
        prefix = op.join(self.args.out_dir, name)
        pat_path = prefix + PAT_SUFF
        unq_path = prefix + UNQ_SUFF
        os.system('cat ' + ' '.join([p for p, u in res]) + ' > ' + pat_path)  # pat
        if self.args.unq:
            os.system('cat ' + ' '.join([u for p, u in res]) + ' > ' + unq_path)  # unq

        # remove all small files
        mult_safe_remove([x[0] for x in res])
        if self.args.unq:
            mult_safe_remove([x[1] for x in res])

        # generate beta file and bgzip the pat, unq files:
        eprint('[bam2pat] bgzipping and indexing:')
        for f in (pat_path, unq_path):
            if not op.isfile(f):
                continue
            subprocess.call(f'bgzip -f@ 14 {f} && tabix -fCb 2 -e 2 {f}.gz', shell=True)
            eprint(f'[bam2pat] generated {f}.gz')

        beta_path = pat2beta(f'{pat_path}.gz', self.out_dir, args=self.args)
        eprint(f'[bam2pat] generated {beta_path}')


def parse_bam2pat_args(parser):
    parser.add_argument('-l', '--lbeta', action='store_true', help='Use lbeta file (uint16) instead of beta (uint8)')
    parser.add_argument('--unq', action='store_true', help='generage unq format as well')


def add_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('bam_path')
    add_GR_args(parser)
    parser.add_argument('--out_dir', '-o', default='.')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--blueprint', '-bp', action='store_true',
            help='filter bad BS conversion reads if <90 percent of CHs are converted')
    parser.add_argument('-F', '--exclude_flags', type=int,
                        help='flags to exclude from bam file (samtools view parameter) ' \
                             '[{}]'.format(FLAGS_FILTER), default=FLAGS_FILTER)
    parser.add_argument('-q', '--mapq', type=int,
                        help='Minimal mapping quality (samtools view parameter) [{}]'.format(MAPQ),
                        default=MAPQ)
    parser.add_argument('-T', '--temp_dir', help='passed to unix sort. Useful in case bam file is very large')
    lists = parser.add_mutually_exclusive_group()
    lists.add_argument('--blacklist', help='bed file. Ignore reads overlapping this bed file',
                        nargs='?', const=True, default=False)
    lists.add_argument('-L', '--whitelist', help='bed file. Consider only reads overlapping this bed file',
                        nargs='?', const=True, default=False)
    add_multi_thread_args(parser)

    return parser


def parse_args(parser):
    parse_bam2pat_args(parser)
    args = parser.parse_args()
    return args


def main():
    """
    Run the WGBS pipeline to generate pat, unq, beta files out of an input bam file
    """
    parser = add_args()
    args = parse_args(parser)
    Bam2Pat(args)


if __name__ == '__main__':
    main()
