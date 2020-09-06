#!/usr/bin/python3 -u


import os
import os.path as op
import argparse
import subprocess
import re
from multiprocessing import Pool
from utils_wgbs import IllegalArgumentError, match_maker_tool, patter_tool, add_GR_args, eprint, add_multi_thread_args
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
    else:
        os.system(cmd)
        return
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()
        if p.returncode or not output:
            eprint(cmd)
            eprint("Failed with subprocess %d\n%s\n%s" % (p.returncode, output.decode(), error.decode()))
            raise IllegalArgumentError('Failed')


def pat_unq(out_path, debug):
    try:
        # sort
        tmp_path = out_path + '.tmp'

        cmd = "sort " + out_path + " -k2,2n -k3,3 -o " + tmp_path
        subprocess_wrap(cmd, debug)

        # break output file into pat and unq:
        # pat file:
        pat_path = out_path + PAT_SUFF
        cmd = 'awk \'{print $1,$2,$3}\' ' + tmp_path + ' | uniq -c | awk \'{OFS="\\t"; print $2,$3,$4,$1}\' > ' + pat_path
        subprocess_wrap(cmd, debug)

        # unq file:
        unq_path = out_path + UNQ_SUFF
        cmd = "sort {} -k4,4n -k3,3 -o {}".format(out_path, tmp_path)
        subprocess_wrap(cmd, debug)
        cmd = 'awk \'{print $1,$4,$5,$3}\' ' + tmp_path + ' | uniq -c | awk \'{OFS="\\t"; print $2,$3,$4,$5,$1}\' > ' + \
              unq_path
        subprocess_wrap(cmd, debug)

        if not debug:
            os.remove(out_path)
            os.remove(tmp_path)

        return pat_path, unq_path
    except IllegalArgumentError as e:
        return None


def proc_chr(input_path, out_path, region, genome, paired_end, ex_flags, mapq, debug):
    """ Convert a temp single chromosome file, extracted from a bam file,
        into two output files: pat and unq."""

    # Run patter tool on a single chromosome. out_path will have the following fields:
    # chr   CpG   Pattern   begin_loc   length(bp)

    # use samtools to extract only the reads from 'chrom'
    flag = '-f 3' if paired_end else ''
    cmd = "samtools view {} {} -q {} -F {} {} | ".format(input_path, region, mapq, ex_flags, flag)
    if debug:
        cmd += ' head -200 | '
    if paired_end:
        # change reads order, s.t paired reads will appear in adjacent lines
        cmd += "{} | ".format(match_maker_tool)
    cmd += "{} {} {} > {}".format(patter_tool, genome.genome_path, genome.chrom_cpg_sizes, out_path)
    # print(cmd)
    subprocess_wrap(cmd, debug)

    return pat_unq(out_path, debug)


class Bam2Pat:
    def __init__(self, args):
        self.args = args
        self.out_dir = args.out_dir
        self.bam_path = args.bam_path
        self.debug = args.debug
        self.gr = GenomicRegion(args)
        self.validate_input()

    def validate_input(self):

        # validate bam path:
        eprint('bam:', self.bam_path)
        if not (op.isfile(self.bam_path) and self.bam_path.endswith('.bam')):
            raise IllegalArgumentError('Invalid bam: {}'.format(self.bam_path))

        # check if bam is sorted by coordinate:
        peek_cmd = 'samtools view -H {} | head -1'.format(self.bam_path)
        if 'coordinate' not in subprocess.check_output(peek_cmd, shell=True).decode():
            raise IllegalArgumentError('bam file must be sorted by coordinate')

        # check if bam is indexed:
        if not (op.isfile(self.bam_path + '.bai')):
            eprint('bai file was not found! Generating...')
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
            cmd = 'samtools idxstats {} | cut -f1 '.format(self.bam_path)
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = p.communicate()
            if p.returncode or not output:
                eprint(cmd)
                eprint("Failed with samtools idxstats %d\n%s\n%s" % (p.returncode, output.decode(), error.decode()))
                eprint('falied to find chromosomes')
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

    def start_threads(self):
        """ Parse each chromosome file in a different process,
            and concatenate outputs to pat and unq files """

        name = op.join(self.out_dir, op.basename(self.bam_path)[:-4])
        processes = []
        with Pool(self.args.threads) as p:
            for c in self.set_regions():
                out_path = name + '_' + c + '.output.tmp'
                params = (self.bam_path, out_path, c, self.gr.genome,
                          self.is_pair_end(), self.args.exclude_flags,
                          self.args.mapq, self.debug)
                processes.append(p.apply_async(proc_chr, params))
            if not processes:
                raise IllegalArgumentError('Empty bam file')
            p.close()
            p.join()
        res = [pr.get() for pr in processes]  # [(pat_path, unq_path) for each chromosome]
        if None in res:
            eprint('threads failed')
            return

        # Concatenate chromosome files
        pat_path = name + PAT_SUFF
        unq_path = name + UNQ_SUFF
        os.system('cat ' + ' '.join([p for p, u in res]) + ' > ' + pat_path)  # pat
        os.system('cat ' + ' '.join([u for p, u in res]) + ' > ' + unq_path)  # unq

        # remove all small files
        list(map(os.remove, [x for l in res for x in l]))

        # generate beta file and bgzip the pat, unq files:
        beta_path = pat2beta(pat_path, self.out_dir, args=self.args)
        eprint('bgzipping and indexing:')
        for f in (pat_path, unq_path):
            eprint('{}...'.format(f))
            subprocess.call('bgzip -f@ 14 {f} && tabix -fCb 2 -e 2 {f}.gz'.format(f=f), shell=True)


def parse_bam2pat_args(parser):
    parser.add_argument('-l', '--lbeta', action='store_true', help='Use lbeta file (uint16) instead of beta (uint8)')


def add_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('bam_path')
    add_GR_args(parser)
    parser.add_argument('--out_dir', '-o', default='.')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('-F', '--exclude_flags', type=int,
                        help='flags to exclude from bam file (samtools view parameter) ' \
                             '[{}]'.format(FLAGS_FILTER), default=FLAGS_FILTER)
    parser.add_argument('-q', '--mapq', type=int,
                        help='Minimal mapping quality (samtools view parameter) [{}]'.format(MAPQ),
                        default=MAPQ)
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
    Bam2Pat(args).start_threads()


if __name__ == '__main__':
    main()
