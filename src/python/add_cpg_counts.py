#!/usr/bin/python3 -u

import os
import os.path as op
import subprocess
import shlex
import datetime
from multiprocessing import Pool
import argparse
from init_genome import chromosome_order
from utils_wgbs import IllegalArgumentError, match_maker_tool, eprint, \
        add_cpg_count_tool, validate_local_exe, add_GR_args, add_multi_thread_args, \
        safe_remove, check_samtools_version, pretty_name
from bam2pat import subprocess_wrap, validate_bam, is_pair_end, MAPQ, \
        FLAGS_FILTER, add_samtools_view_flags, is_region_empty, set_regions
from genomic_region import GenomicRegion
from convert import load_bed
import pandas as pd
import numpy as np
import os.path as op
import uuid


BAM_SUFF = '.bam'

# Minimal Mapping Quality to consider.
# 10 means include only reads w.p. >= 0.9 to be mapped correctly.
# And missing values (255)


def proc_chr(input_path, out_path_name, region, genome, paired_end, ex_flags, mapq,
             debug, verbose, min_cpg, clip, bed_file, extended_bed, add_pat, in_flags,
             drop_singles):
    """ Convert a temp single chromosome file, extracted from a bam file,
        into a sam formatted (no header) output file."""

    # Run patter tool 'bam' mode on a single chromosome

    unsorted_bam = out_path_name + '_unsorted.output.bam'
    out_path = out_path_name + '.output.bam'
    out_directory = os.path.dirname(out_path)

    # use samtools to extract only the reads from 'chrom'
    # flag = '-f 3' if paired_end else ''
    if in_flags is None:
        in_flags = '-f 3' if paired_end else ''
    else:
        in_flags = f'-f {in_flags}'
    cmd = f'samtools view {input_path} {region} -q {mapq} -F {ex_flags} {in_flags} -P'
    if bed_file is not None:
        cmd += f' -L {bed_file} '
    # check if chrom is empty:
    if is_region_empty(cmd, region, verbose=verbose):
        return ''

    if debug:
        cmd += ' | head -200 '
    if paired_end:
        # change reads order, s.t paired reads will appear in adjacent lines
        ds = ' --drop_singles' if drop_singles else ''
        cmd += f' | {match_maker_tool}{ds} '
    cmd += f' | {add_cpg_count_tool} {genome.dict_path} {region} --clip {clip}'
    if bed_file is not None:
        cmd += f' --bed_file {bed_file} --bed_ext_file {extended_bed}'
    if add_pat:
        cmd += ' --pat'
    cmd += f' --min_cpg {min_cpg}'
    cmd += f' | cat <( samtools view -H {input_path} ) - | samtools view -hb - > {unsorted_bam}'

    sort_cmd = f'samtools sort -o {out_path} -T {out_directory} {unsorted_bam}'  # TODO: use temp directory, as in bam2pat

    subprocess_wrap(cmd, verbose)
    subprocess_wrap(sort_cmd, verbose)
    safe_remove(unsorted_bam)
    return out_path


class BamMethylData:
    def __init__(self, args, bam_path):
        self.args = args
        self.out_dir = args.out_dir
        self.bam_path = bam_path
        self.bed_path = args.bed_file
        self.gr = GenomicRegion(args)
        self.extended_bed_path = self.validate_bed()
        self.validate_input()

    def validate_bed(self):
        if self.bed_path is None:
            return
        df = load_bed(self.bed_path).iloc[:, :3]

        # check start before end
        if ((df['start'] >= df['end']).any()):
            eprint(f'[wt add_cpg_counts] bed file is not legal - end before start')
            raise IllegalArgumentError('Bed file is not legal')

        # check start is monotonic
        ref_chroms =  self.gr.genome.get_chroms()
        for chrom in ref_chroms:
            if not df[df['chr'] == chrom]['start'].is_monotonic_increasing:
                eprint(f'[wt add_cpg_counts] bed file is not sorted')
                raise IllegalArgumentError('Bed file is not sorted')

        # create tmp files from bed file (extend each region with +-1000 bp)
        df['start'] = np.maximum(1, df['start'] - 1000)
        df['end'] += 1000
        extended_bed_path = op.join(self.out_dir, pretty_name(self.bam_path) + str(uuid.uuid4()))[:6] + '.tmp.bed'
        df.to_csv(extended_bed_path, sep='\t', header=False, index=False)
        return extended_bed_path

    def validate_input(self):
        # validate output dir:
        if not (op.isdir(self.out_dir)):
            raise IllegalArgumentError(f'Invalid output dir: {self.out_dir}')

    def start_threads(self):
        """ Parse each chromosome file in a different process,
            and concatenate outputs to pat and unq files """
        print(datetime.datetime.now().isoformat() + ': *** starting processing of each chromosome')
        name = op.join(self.out_dir, pretty_name(self.bam_path))
        is_PE = is_pair_end(self.bam_path, self.gr.genome)
        if self.gr.region_str is None:
            final_path = name + f'.{self.args.suffix}' + BAM_SUFF
            params = []
            for c in set_regions(self.bam_path, self.gr):
                out_path_name = name + '_' + c
                params.append((self.bam_path, out_path_name, c, self.gr.genome,
                        is_PE, self.args.exclude_flags, self.args.mapq,
                        self.args.debug, self.args.verbose, self.args.min_cpg,
                        self.args.clip, self.args.bed_file, self.extended_bed_path,
                        self.args.add_pat, self.args.include_flags, self.args.drop_singles))
            p = Pool(self.args.threads)
            res = p.starmap(proc_chr, params)
            p.close()
            p.join()
        else:
            region_str_for_name = self.gr.region_str.replace(':', '_').replace('-', '_')
            final_path = name + f'.{region_str_for_name}.{self.args.suffix}{BAM_SUFF}'
            out_path_name = name + '_' + '1'
            res = [proc_chr(self.bam_path, out_path_name, self.gr.region_str, self.gr.genome,
                            is_PE, self.args.exclude_flags, self.args.mapq,
                            self.args.debug, self.args.verbose, self.args.min_cpg, self.args.clip,
                            self.args.bed_file, self.extended_bed_path, self.args.add_pat,
                            self.args.include_flags, self.args.drop_singles)]
        print('finished adding CpG counts')
        if None in res:
            print('threads failed')
            return

        res = [r for r in res.copy() if r != '']
        print(datetime.datetime.now().isoformat() + ': finished processing each chromosome')
        if not res:
            eprint(f'[wt add_cpg_counts] no reads found for {self.bam_path}.')
            eprint(f'                    run wgbstools add_cpg_counts with --verbose for more information')
            return
        # Concatenate chromosome files

        out_directory = os.path.dirname(final_path)
        cmd = f'samtools merge -c -p -f {final_path} ' + ' '.join([p for p in res])
        print(datetime.datetime.now().isoformat() + ': starting cat of files')
        process = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        stdout, stderr = process.communicate()
        print(datetime.datetime.now().isoformat() + ': finished cat of files')

        idx_command = f'samtools index {final_path}'
        print('starting index of output bam ' + datetime.datetime.now().isoformat())
        idx_process = subprocess.Popen(shlex.split(idx_command), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        stdout, stderr = idx_process.communicate()
        print(datetime.datetime.now().isoformat() + ': finished index of output bam')
        # remove all small files
        list(map(os.remove, [l for l in res]))
        # remove tmp files
        safe_remove(self.extended_bed_path) # the file that we create from bed file

def add_cpg_args(parser):
    parser.add_argument('--drop_singles',  action='store_true',
                        help='For paired bam only - if mate is not exists drop single read')
    parser.add_argument('--suffix', default='counts',
                        help='The output file suffix. The output file will be [in_file].[suffix].bam. By default the '
                             'suffix is "counts".')
    parser.add_argument('--add_pat', action='store_true',
                        help='Indicates whether to add the methylation pattern of the read (pair).')
    return parser

def add_args(parser):
    parser.add_argument('bam', nargs='+')
    add_GR_args(parser, bed_file=True)
    parser.add_argument('--out_dir', '-o', default='.')
    parser.add_argument('--min_cpg', type=int, default=1,
                help='Reads covering less than MIN_CPG sites are removed [1]')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--force', '-f', action='store_true', help='overwrite existing files if exists')
    parser.add_argument('--verbose', '-v', action='store_true')
    add_samtools_view_flags(parser)
    parser.add_argument('--clip', type=int, default=0,
                        help='Clip for each read the first and last CLIP characters [0]')
    add_multi_thread_args(parser)

    return parser


def main():
    """
    Add to bam file an extra field, YI:Z:{nr_meth},{nr_unmeth},
    to count Cytosine retention at CpG context.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser = add_args(parser)
    parser = add_cpg_args(parser)
    args = parser.parse_args()
    validate_local_exe(add_cpg_count_tool)
    if not check_samtools_version(minor=15, verbose=args.verbose):
        eprint('[wt add_cpg_counts] Error: add_cpg_counts only works with samtools version >=1.15.\n'
                '                    Please update samtools.\n'
                '                    Run wgbstools add_cpg_counts with --verbose for more information')
        return
    for bam in args.bam:
        if not validate_bam(bam):
            eprint(f'[wt add_cpg_counts] Skipping {bam}')
            continue
        BamMethylData(args, bam).start_threads()


if __name__ == '__main__':
    main()
