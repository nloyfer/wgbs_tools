#!/usr/bin/python3 -u

import os
import os.path as op
import argparse
import subprocess
from multiprocessing import Pool
from utils_wgbs import IllegalArgumentError, match_maker_tool, patter_tool, add_GR_args
from init_genome_ref_wgbs import chromosome_order
from pat2beta import pat2beta
from pipeline_wgbs.test import run_test
from genomic_region import GenomicRegion


PAT_SUFF = '.pat'
UNQ_SUFF = '.unq'

# Minimal Mapping Quality to consider.
# 10 means include only reads w.p. >= 0.9 to be mapped correctly.
# And missing values (255)
MAPQ = 10

FLAGS_FILTER = 1796  # filter flags with these bits
# todo: unsorted / sorted by name


def pat_unq(out_path):
    # sort and drop duplicated lines # todo: need to drop duplicates?
    tmp_path = out_path + '.tmp'

    subprocess.call("sort " + out_path + " -k2,2n -k3,3 -o " + tmp_path, shell=True)

    # break output file into pat and unq:
    # pat file:
    pat_path = out_path + PAT_SUFF
    cmd = 'awk \'{print $1,$2,$3}\' ' + tmp_path + ' | uniq -c | awk \'{OFS="\\t"; print $2,$3,$4,$1}\' > ' + pat_path
    subprocess.call(cmd, shell=True)

    # unq file:
    unq_path = out_path + UNQ_SUFF
    subprocess.call("sort {} -k4,4n -k3,3 -o {}".format(out_path, tmp_path), shell=True)
    cmd = 'awk \'{print $1,$4,$5,$3}\' ' + tmp_path + ' | uniq -c | awk \'{OFS="\\t"; print $2,$3,$4,$5,$1}\' > ' +\
          unq_path
    subprocess.call(cmd, shell=True)

    os.remove(out_path)
    os.remove(tmp_path)
    return pat_path, unq_path


def proc_chr(input_path, out_path, region, genome, paired_end, debug):
    """ Convert a temp single chromosome file, extracted from a bam file,
        into two output files: pat and unq."""

    # Run patter tool on a single chromosome. out_path will have the following fields:
    # chr   CpG   Pattern   begin_loc   length(bp)

    # use samtools to extract only the reads from 'chrom'
    flag = '-f 3' if paired_end else ''
    cmd = "samtools view {} {} -q {} -F 1796 {} | ".format(input_path, region, MAPQ, flag)
    if debug:
        cmd += ' head -200 | '
    if paired_end:
        # change reads order, s.t paired reads will appear in adjacent lines
        cmd += "{} | ".format(match_maker_tool)
    cmd += "{} {} {} > {}".format(patter_tool, genome.genome_path, genome.chrom_cpg_sizes, out_path)  # todo: pipe to sort? check if faster
    print(cmd)
    subprocess.call(cmd, shell=True)

    return pat_unq(out_path)


class Bam2Pat:
    def __init__(self, bam_path, out_dir, gr, debug):
        self.out_dir = out_dir
        self.bam_path = bam_path
        self.debug = debug
        self.gr = gr
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

    def is_pair_end(self):
        first_line = subprocess.check_output('samtools view {} | head -1'.format(self.bam_path), shell=True)
        return int(first_line.decode().split('\t')[1]) & 1

    def set_regions(self):
        if self.gr.region_str:
            return [self.gr.region_str]
        else:
            cmd = 'samtools idxstats {} | cut -f1 | grep chr'.format(self.bam_path)
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = p.communicate()
            if p.returncode or not output:
                print("Failed with samtools idxstats %d\n%s\n%s" % (p.returncode, output, error))
                return []
            chroms = list(sorted(output.decode()[:-1].split('\n'), key=chromosome_order))
            return chroms
            # return ['chr{}'.format(i) for i in list(range(1, 23)) + ['X', 'Y', 'M']]

    def start_threads(self):
        """ Parse each chromosome file in a different process,
            and concatenate outputs to pat and unq files """

        name = op.join(self.out_dir, op.basename(self.bam_path)[:-4])
        processes = []
        with Pool() as p:
            for c in self.set_regions():
                out_path = name + '_' + c + '.output.tmp'
                params = (self.bam_path, out_path, c, self.gr.genome, self.is_pair_end(), self.debug)
                processes.append(p.apply_async(proc_chr, params))
            if not processes:
                raise IllegalArgumentError('Empty bam file')
            p.close()
            p.join()
        res = [pr.get() for pr in processes]    # [(pat_path, unq_path) for each chromosome]

        # Concatenate chromosome files
        pat_path = name + PAT_SUFF
        unq_path = name + UNQ_SUFF
        os.system('cat ' + ' '.join([p for p, u in res]) + ' > ' + pat_path)  # pat
        os.system('cat ' + ' '.join([u for p, u in res]) + ' > ' + unq_path)  # unq

        # remove all small files
        list(map(os.remove, [x for l in res for x in l]))

        # generate beta file and bgzip the pat, unq files:
        beta_path = pat2beta(pat_path, self.out_dir)
        print('bgzipping...')
        for f in (pat_path, unq_path):
            subprocess.call('bgzip -f@ 14 {f} && tabix -fCb 2 -e 2 {f}.gz'.format(f=f), shell=True)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('bam_path')
    add_GR_args(parser)
    parser.add_argument('--out_dir', '-o', default='.')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--test', action='store_true',
                        help='Perform a test for the pipeline. Ignore other parameters.')
    args = parser.parse_args()
    return args


def main():
    """
    Run the WGBS pipeline to generate pat, unq, beta files out of an input bam file
    """
    args = parse_args()

    if args.test:
        print('Testing...')
        run_test.main()
        return

    # else
    Bam2Pat(args.bam_path, args.out_dir, GenomicRegion(args), args.debug).start_threads()


if __name__ == '__main__':
    main()
