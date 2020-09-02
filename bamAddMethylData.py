#!/usr/bin/python3 -u

import os
import os.path as op
import argparse
import subprocess
import shlex
import re
import datetime
import multiprocessing
from multiprocessing import Pool
from utils_wgbs import IllegalArgumentError, patter_tool, match_maker_tool, add_GR_args, eprint
from bam2pat import parse_args, subprocess_wrap, CHROMS
from init_genome_ref_wgbs import chromosome_order
from genomic_region import GenomicRegion


BAM_SUFF = '_md.bam'

# Minimal Mapping Quality to consider.
# 10 means include only reads w.p. >= 0.9 to be mapped correctly.
# And missing values (255)
MAPQ = 10

# todo: unsorted / sorted by name


def proc_chr(input_path, out_path_name, region, genome, header_path, paired_end, ex_flags, mapq, debug):
    """ Convert a temp single chromosome file, extracted from a bam file,
        into a sam formatted (no header) output file."""

    # Run patter tool 'bam' mode on a single chromosome

    unsorted_bam = out_path_name + '_unsorted.output.bam'
    out_path = out_path_name + '.output.bam'
    out_directory = os.path.dirname(out_path)

    # use samtools to extract only the reads from 'chrom'
    flag = '-f 3'
    cmd = "samtools view {} {} -q {} -F {} {} | ".format(input_path, region, mapq, ex_flags, flag)
    if debug:
        cmd += ' head -200 | '
    if paired_end:
        # change reads order, s.t paired reads will appear in adjacent lines
        cmd += "{} | ".format(match_maker_tool)
    cmd += "{} {} {} bam | cat {} - > {}".format(patter_tool, genome.genome_path,
                                                                            genome.chrom_cpg_sizes,
                                                                            header_path, unsorted_bam)

    sort_cmd = 'samtools sort -o {} -T {} {}'.format(out_path, out_directory, unsorted_bam)

    #print(cmd)
    subprocess_wrap(cmd, debug)
    subprocess_wrap(sort_cmd, debug)
    os.remove(unsorted_bam)
    return out_path

def get_header_command(input_path):
    return "samtools view -H {}".format(input_path)

def proc_header(input_path, out_path, debug):
    """ extracts header from bam file and saves it to tmp file."""

    cmd = get_header_command(input_path) + " > {} ".format(out_path)
    #print(cmd)
    subprocess_wrap(cmd, debug)

    return out_path


class BamMethylData:
    def __init__(self, args):
        self.args = args
        self.out_dir = args.out_dir
        self.bam_path = args.bam_path
        self.debug = args.debug
        self.gr = GenomicRegion(args)
        self.validate_input()

    def validate_input(self):

        # validate bam path:
        print('bam:', self.bam_path)
        if not (op.isfile(self.bam_path) and self.bam_path.endswith('.bam')):
            raise IllegalArgumentError('Invalid bam: {}'.format(self.bam_path))

        # check if bam is sorted by coordinate:
        peek_cmd = 'samtools view -H {} | head -1'.format(self.bam_path)
        if 'coordinate' not in subprocess.check_output(peek_cmd, shell=True).decode():
            raise IllegalArgumentError('bam file must be sorted by coordinate')

        # check if bam is indexed:
        if not (op.isfile(self.bam_path + '.bai')):
            print('bai file was not found! Generating...')
            r = subprocess.call(['samtools', 'index', self.bam_path])
            if r:
                raise IllegalArgumentError('Failed indexing bam: {}'.format(self.bam_path))

        # validate output dir:
        if not (op.isdir(self.out_dir)):
            raise IllegalArgumentError('Invalid output dir: {}'.format(self.out_dir))

    def set_regions(self):
        if self.gr.region_str:
            return [self.gr.region_str]
        else:
            cmd = 'samtools idxstats {} | cut -f1 '.format(self.bam_path)
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = p.communicate()
            if p.returncode or not output:
                print(cmd)
                print("Failed with samtools idxstats %d\n%s\n%s" % (p.returncode, output.decode(), error.decode()))
                print('falied to find chromosomes')
                return []
            nofilt_chroms = output.decode()[:-1].split('\n')
            filt_chroms = [c for c in nofilt_chroms if 'chr' in c]
            if not filt_chroms:
                filt_chroms = [c for c in nofilt_chroms if c in CHROMS]
            else:
                filt_chroms = [c for c in filt_chroms if re.match(r'^chr([\d]+|[XYM])$', c)]
            chroms = list(sorted(filt_chroms, key=chromosome_order))
            if not chroms:
                eprint('Failed retrieving valid chromosome names')
                raise IllegalArgumentError('Failed')

            return chroms

    def is_pair_end(self):
        first_line = subprocess.check_output('samtools view {} | head -1'.format(self.bam_path), shell=True)
        return int(first_line.decode().split('\t')[1]) & 1

    def intermediate_bam_file_view(self, name):
        return '<(samtools view {})'.format(name)

    def process_substitute(self, cmd):
        return '<({})'.format(cmd)

    def start_threads(self):
        """ Parse each chromosome file in a different process,
            and concatenate outputs to pat and unq files """
        print(datetime.datetime.now().isoformat() + ": *** starting processing of each chromosome")
        name = op.join(self.out_dir, op.basename(self.bam_path)[:-4])
        header_path = name + '.header'
        proc_header(self.bam_path, header_path, self.debug)
        processes = []
        with Pool(self.args.threads) as p:
            for c in self.set_regions():
                out_path_name = name + '_' + c
                params = (self.bam_path, out_path_name, c, self.gr.genome,
                        header_path, self.is_pair_end(), self.args.exclude_flags,
                        self.args.mapq, self.debug)
                processes.append(p.apply_async(proc_chr, params))
            if not processes:
                raise IllegalArgumentError('Empty bam file')
            p.close()
            p.join()
        res = [pr.get() for pr in processes]    # [(pat_path, unq_path) for each chromosome]
        print('finished patter')
        if None in res:
            print('threads failed')
            return

        print(datetime.datetime.now().isoformat() + ": finished processing each chromosome")
        # Concatenate chromosome files
        final_path = name + BAM_SUFF
        out_directory = os.path.dirname(final_path)
        # cmd = '/bin/bash -c "cat <({})'.format(get_header_command(self.bam_path)) + ' ' +\
        #       ' '.join([self.intermediate_bam_file_view(p) for p in res]) + ' | samtools view -b - > ' + final_path_unsorted + '"'
        cmd ="samtools merge -h {} {} ".format(header_path, final_path) + ' '.join([p for p in res])
        # cmd = '/bin/bash -c "samtools cat -h <({})'.format(get_header_command(self.bam_path)) + ' ' + \
        #       ' '.join(
        #           [p for p in res]) + ' > ' + final_path + '"'
        print(datetime.datetime.now().isoformat() + ': starting cat of files')
        process = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        stdout, stderr = process.communicate()
        print(datetime.datetime.now().isoformat() + ": finished cat of files")

        # sort_cmd = 'samtools sort -o {} -T {} {}'.format(final_path, out_directory, final_path_unsorted)
        # print(datetime.datetime.now().isoformat() + ': starting sort of file')
        # sort_process = subprocess.Popen(shlex.split(sort_cmd), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        # stdout, stderr = sort_process.communicate()
        # print(datetime.datetime.now().isoformat() + ": finished sort of file")

        idx_command = "samtools index {}".format(final_path)
        print('starting index of output bam ' + datetime.datetime.now().isoformat())
        idx_process = subprocess.Popen(shlex.split(idx_command), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        stdout, stderr = idx_process.communicate()
        print(datetime.datetime.now().isoformat() + ": finished index of output bam")
        res.append(header_path)
        # remove all small files
        list(map(os.remove, [l for l in res]))


def main():
    """
    Add to bam file an extra field, YI:Z:{nr_meth},{nr_unmeth},
    to count Cytosine retention at CpG context.
    """
    args = parse_args()
    BamMethylData(args).start_threads()


if __name__ == '__main__':
    main()
