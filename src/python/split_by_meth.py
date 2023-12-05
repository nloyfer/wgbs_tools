#!/usr/bin/python3 -u
import os
import os.path as op
import argparse
import subprocess
from multiprocessing import Pool

from genomic_region import GenomicRegion
from utils_wgbs import IllegalArgumentError, eprint, add_multi_thread_args, bam_meth_split_tool, add_GR_args


def subprocess_wrap(cmd, debug):
    if debug:
        print(cmd)
        return
    os.system(cmd)


def proc_chr(bam, name, region, is_meth, homog_prop, min_cpg, ex_flags, mapq, debug, verbose):
    # use samtools to extract only the reads from 'chrom'
    # flag = '-f 3' if paired_end else ''
    prop_to_use = float(homog_prop) if is_meth else (1 - float(homog_prop))
    view_cmd = f"samtools view {bam} "
    if region is not None:
        view_cmd = view_cmd + f"{region} "
    if mapq is not None:
        view_cmd = view_cmd + f"-q {mapq} "
    if ex_flags is not None:
        view_cmd = view_cmd + f"-F {ex_flags} "
    split_cmd = f"{view_cmd} | {bam_meth_split_tool} {prop_to_use} {min_cpg}"
    samtools_header_cmd = f"samtools view -H {bam}"
    cmd = f'/bin/bash -c "cat <({samtools_header_cmd}) <({split_cmd}) | samtools view -h -b - > {name} && samtools' \
          f' index {name}"'

    if debug:
        cmd += ' | head -200 '
    if verbose:
        print(cmd)
    subprocess_wrap(cmd, debug)


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

    # check if bam has YI:Z flag:
    peek_cmd = f'samtools view {bam} | head -1'
    if 'YI:Z' not in subprocess.check_output(peek_cmd, shell=True).decode():
        eprint('bam file must contain CpG counts info. Please run `wgbstools add_cpg_counts`')
        return False

    # check if bam is indexed:
    if not (op.isfile(bam + '.bai')):
        eprint('[wt bam2pat] bai file was not found! Generating...')
        if subprocess.call(['samtools', 'index', bam]):
            eprint(f'[wt bam2pat] Failed indexing bam: {bam}')
            return False
    return True


class MethSplit:
    def __init__(self, args, bam):
        self.args = args
        self.tmp_dir = None
        self.verbose = args.verbose
        self.out_dir = args.out_dir
        self.bam_path = bam
        self.homog_prop = args.homog_prop
        self.min_cpg = args.min_cpg
        self.gr = GenomicRegion(args)
        self.start_threads()

    def start_threads(self):
        base_name = op.basename(self.bam_path)[:-4]
        if self.gr.region_str is not None:
            new_region_str = self.gr.region_str.replace(":", "_").replace("-", "_")
            base_name = base_name + f".{new_region_str}"
        meth_name = op.join(self.out_dir, base_name + ".M.bam")
        unmeth_name = op.join(self.out_dir, base_name + ".U.bam")

        params = [
            (self.bam_path, meth_name, self.gr.region_str, True, self.homog_prop, self.min_cpg, self.args.exclude_flags,
             self.args.mapq, self.args.debug, self.verbose),
            (self.bam_path, unmeth_name, self.gr.region_str, False, self.homog_prop, self.min_cpg,
             self.args.exclude_flags,
             self.args.mapq, self.args.debug, self.verbose)]

        p = Pool(self.args.threads)
        p.starmap(proc_chr, params)

        p.close()
        p.join()


def add_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('bam', help="The full path of the bam file to process", nargs='+')
    parser.add_argument('homog_prop', help="A fraction with which to determine homogenous reads. All reads with "
                                           "methylation proportion >= [homog_prop] will be classified as highly "
                                           "methylated reads while all reads with methylation proportion <= 1 - [homog_prop] "
                                           "will be classified as mostly un-methylated reads. Must be above 0.5")
    parser.add_argument('--min_cpg',
                        help="The value of CpGs required per fragment (read pair combined number of CpGs). "
                             "default [1]",
                        default=1)
    parser.add_argument('--out_dir', '-o', default='.')
    parser.add_argument('--force', '-f', action='store_true', help='overwrite existing files if exists')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('-F', '--exclude_flags', type=int,
                        help='flags to exclude from bam file (samtools view parameter)',
                        default=None)
    parser.add_argument('-q', '--mapq', type=int,
                        help=f'Minimal mapping quality (samtools view parameter)',
                        default=None)
    add_GR_args(parser)
    add_multi_thread_args(parser)

    return parser


def main():
    """
    Split the input bam file according to methylation fraction. That is the input is a bam file (which is the output of
    wgbstools add_cpg_counts), a fraction `homog_prop`, and the output
    is two bam files: one containing reads with methylation proportion above the `homog_prop` and one containing reads
    with methylation proportion below 1-`homog_prop`.
    """
    parser = add_args()
    args = parser.parse_args()
    # validate output dir:
    if not op.isdir(args.out_dir):
        raise IllegalArgumentError(f'Invalid output dir: {args.out_dir}')

    if (float(args.homog_prop) <= 0.5):
        raise IllegalArgumentError(f"Illegal homog_prop value {args.homog_prop}. Must be above 0.5")

    for bam in args.bam:
        if not validate_bam(bam):
            eprint(f'[wt bam2pat] Skipping {bam}')
            continue

        MethSplit(args, bam)


if __name__ == '__main__':
    main()
