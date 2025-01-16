#!/usr/bin/python3 -u


import os
import os.path as op
import argparse
import subprocess
import shutil
import re
from multiprocessing import Pool

from bam2pat import Bam2Pat, add_args, parse_bam2pat_args, extend_region, validate_bam
from utils_wgbs import IllegalArgumentError, match_maker_tool, eprint, \
    add_multi_thread_args, mult_safe_remove, GenomeRefPaths, validate_single_file, \
    delete_or_skip, allele_split_tool, add_no_beta_arg, validate_local_exe, add_no_pat_arg
from init_genome import chromosome_order

PAT_SUFF = '.pat.gz'

# Minimal Mapping Quality to consider.
# 10 means include only reads w.p. >= 0.9 to be mapped correctly.
# And missing values (255)
MAPQ = 10
SNP_Q = 0
FLAGS_FILTER = 1796  # filter flags with these bits

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

        mult_safe_remove([out_path])

        return pat_path
    except IllegalArgumentError as e:
        return None


def proc_chr(bam, out_path, name, snp_pos, snp_let1, snp_let2, ex_flags, mapq, debug, verbose, no_beta, no_pat, snp_qual, genome):
    """ Convert a temp single chromosome file, extracted from a bam file, into pat """

    # Run patter tool on a single chromosome. out_path will have the following fields:
    # chr   CpG   Pattern   begin_loc   length(bp)


    # use samtools to extract only the reads from 'chrom'
    # flag = '-f 3' if paired_end else ''
    chrom, position = snp_pos.split(":")
    region = f"{snp_pos}-{int(position) + 1}"
    region = extend_region(region)
    cmd = f'samtools view {bam} {region} -q {mapq} -F {ex_flags} ' #{flag} '
    if debug:
        cmd += ' | head -200 '
    cmd += f' | {match_maker_tool} - |'
    bam_file_out = op.join(out_path, name) + ".bam"
    cmd += f' {allele_split_tool} --snp_pos {position} --snp_let1 \'{snp_let1}\' --snp_let2 \'{snp_let2}\' ' \
           f'--qual_filter {snp_qual}'
    final_cmd = f'/bin/bash -c "cat <(samtools view -H {bam}) <({cmd}) | samtools view -h -b - | samtools sort -O bam -' \
                f' > {bam_file_out} && samtools index {bam_file_out}"'
    if verbose:
        print(cmd)
    subprocess_wrap(final_cmd, debug)
    if not no_pat:
        parser = argparse.ArgumentParser()
        parser = add_args(parser)
        parse_bam2pat_args(parser)
        bam2patargs_list = [bam_file_out, "--out_dir", out_path, "-r", chrom, "--threads", "1", '--genome', genome]
        if no_beta:
            bam2patargs_list += ["--no_beta"]
        args = parser.parse_args(bam2patargs_list)
        Bam2Pat(args, bam_file_out)


def is_pair_end(bam):
    first_line = subprocess.check_output(f'samtools view {bam} | head -1', shell=True)
    return int(first_line.decode().split('\t')[1]) & 1


class SNPSplit:
    def __init__(self, args, bam):
        self.args = args
        self.tmp_dir = None
        self.verbose = args.verbose
        self.out_dir = args.out_dir
        self.bam_path = bam
        self.snp_pos = args.pos
        self.snp_let1 = args.alleles.split("/")[0]
        self.snp_let2 = args.alleles.split("/")[1]
        self.no_beta = args.no_beta
        self.no_pat = args.no_pat
        self.snp_qual = args.snp_qual
        self.genome = args.genome
        # self.gr = GenomicRegion(args)
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
        name1 = op.basename(self.bam_path)[:-4] + f".{self.snp_pos}" + f".{self.snp_let1}"
        name2 = op.basename(self.bam_path)[:-4] + f".{self.snp_pos}" + f".{self.snp_let2}"

        params = [(self.bam_path, self.out_dir, name1, self.snp_pos, self.snp_let1, self.snp_let2, self.args.exclude_flags,
                   self.args.mapq, self.args.debug, self.verbose, self.no_beta, self.no_pat, self.snp_qual, self.genome),
                  (self.bam_path, self.out_dir, name2, self.snp_pos,
                   self.snp_let2, self.snp_let1,
                   self.args.exclude_flags,
                   self.args.mapq, self.args.debug, self.verbose, self.no_beta, self.no_pat, self.snp_qual, self.genome)]

        p = Pool(self.args.threads)
        p.starmap(proc_chr, params)
        p.close()
        p.join()
        #
        # self.concat_parts(name, pat_parts)
        # # remove all small files
        # mult_safe_remove(pat_parts)
        x = 0


def add_args_snp_splitt():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('bam', help="The full path of the bam file to process")
    parser.add_argument('pos', help="The position of the nucleotide from which to split the bam file. The format of "
                                    "this argument should be 'chrXX:YYYY' where XX is the chromosome identifier and "
                                    "YYYY is the position within the chromosome.")
    parser.add_argument('alleles', help="The genoypes of the alleles on which to split the bam file in the format 'X/Y' where X and Y are characters 'A', 'C', 'G', or 'T' ")
    add_no_beta_arg(parser)
    add_no_pat_arg(parser)
    parser.add_argument('--out_dir', '-o', default='.')
    parser.add_argument('--force', '-f', action='store_true', help='overwrite existing files if exists')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('-F', '--exclude_flags', type=int,
                        help='flags to exclude from bam file (samtools view parameter) ' \
                             f'[{FLAGS_FILTER}]', default=FLAGS_FILTER)
    parser.add_argument('-q', '--mapq', type=int,
                        help=f'Minimal mapping quality (samtools view parameter) [{MAPQ}]',
                        default=MAPQ)
    parser.add_argument('--snp_qual', type=int,
                        help=f'Minimal mapping quality (phred) of base call at the position of the polymorphism [{SNP_Q}]',
                        default=SNP_Q)
    parser.add_argument('--genome', help='Genome reference name. Default is "default".', default='default')
    add_multi_thread_args(parser)

    return parser


def parse_args_snp_split(parser):
    # parse_bam2pat_args(parser)
    args = parser.parse_args()
    return args


def main():
    """
    Run the WGBS pipeline to generate pat & beta files out of an input bam file
    """
    parser = add_args_snp_splitt()
    args = parse_args_snp_split(parser)
    # validate output dir:
    if not op.isdir(args.out_dir):
        raise IllegalArgumentError(f'Invalid output dir: {args.out_dir}')

    if args.no_pat and not args.no_beta:
        print(f'Cannot create a beta file without creating a pat file. If you would like to create a beta file please'
              f' remove the `no_pat` argument.')

    validate_local_exe(allele_split_tool)
    for bam in [args.bam]:
        if not validate_bam(bam):
            eprint(f'[wt bam2pat] Skipping {bam}')
            continue


        snp_let1 = args.alleles.split("/")[0]
        snp_let2 = args.alleles.split("/")[1]
        snp_pos = args.pos
        pat1 = op.join(args.out_dir, op.basename(bam)[:-4] + f".{snp_pos}.{snp_let1}" + PAT_SUFF)
        pat2 = op.join(args.out_dir, op.basename(bam)[:-4] + f".{snp_pos}.{snp_let2}" + PAT_SUFF)
        if not delete_or_skip(pat1, args.force):
            continue
        if not delete_or_skip(pat2, args.force):
            continue
        SNPSplit(args, bam)


if __name__ == '__main__':
    main()
