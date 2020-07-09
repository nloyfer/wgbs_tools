import re
import subprocess
from utils_wgbs import IllegalArgumentError, GenomeRefPaths, eprint
import sys
import numpy as np
import pandas as pd
import os.path as op


class GenomicRegion:
    def __init__(self, args=None, region=None, sites=None, genome_name='hg19'):
        self.genome_name = genome_name
        self.chrom = None
        self.sites = sites
        self.region_str = region
        self.bp_tuple = None
        self.chrs_sz = None  # DataFrame of chromosomes sizes (in number of sites)
        self.args = args

        # todo: this could be prettier
        if args is not None:
            self.genome_name = args.genome
            self.genome = GenomeRefPaths(self.genome_name)
            if args.sites:
                self.parse_sites(args.sites)
            elif args.region:
                self.parse_region(args.region)
        elif region is not None:
            self.genome = GenomeRefPaths(self.genome_name)
            self.parse_region(region)
        elif sites is not None:
            self.genome = GenomeRefPaths(self.genome_name)
            self.parse_sites(sites)
        else:
            raise IllegalArgumentError('Invalid GR init {}'.format(region))

        self.nr_sites = None if self.sites is None else self.sites[1] - self.sites[0]
        self.annotation = self.add_anno()

    def add_anno(self):
        if self.args is None or self.is_whole():
            return
        if self.args.no_anno:
            return
        anno_path = self.genome.annotations
        if anno_path is None:
            return
        try:

            cmd = 'tabix {} {} | cut -f4- | uniq'.format(anno_path, self.region_str)
            res = subprocess.check_output(cmd, shell=True).decode().strip()
            return res
        except subprocess.CalledProcessError:
            eprint('Failed to retrieve annotation for reagion ', self.region_str)
            return

    def parse_sites(self, sites_str):
        """ Parse input of the type -s / --sites (e.g 15-25) """

        # Parse sites string:
        s1, s2 = self._sites_str_to_tuple(sites_str)

        # Translate sites indexes to genomic loci:
        self.chrom, region_from = self.index2locus(s1)
        chrom2, region_to = self.index2locus(s2 - 1)  # non-inclusive
        region_to += 1  # include the whole last site (C and G)
        if self.chrom != chrom2:
            eprint('ERROR: sites range cross chromosomes! ({}, {})'.format(s1, s2))
            raise IllegalArgumentError('Invalid sites input')

        # Update GR fields:
        self.sites = (s1, s2)
        self.region_str = "{}:{}-{}".format(self.chrom, region_from, region_to)
        self.bp_tuple = (region_from, region_to)

    def _chrome_size(self):
        df = pd.read_csv(self.genome.chrom_sizes, sep='\t', header=None, names=['chr', 'size'])
        return int(df[df['chr'] == self.chrom]['size'])

    def parse_region(self, region):
        """ Parse input of the type -r / --region (e.g chr11:200-300) """
        region = region.replace(',', '')  # remove commas
        chrome_match = re.match(r'^chr([\d]+|[XYM])$', region)
        region_match = re.match(r'chr([\d]+|[XYM]):([\d]+)-([\d]+)', region)

        # In case region is a whole chromosome
        if chrome_match:
            self.chrom = 'chr' + chrome_match.group(1)
            region_from = 1
            region_to = self._chrome_size()

        # match region string to format chrom:from-to
        elif region_match:
            self.chrom = 'chr' + region_match.group(1)
            region_from = int(region_match.group(2))
            region_to = int(region_match.group(3))
            if region_to <= region_from:
                raise IllegalArgumentError('Invalid genomic region: {}. end before start'.format(region))
            if region_to > self._chrome_size() or region_from < 1:
                raise IllegalArgumentError('Invalid genomic region: {}. Out of range'.format(region))

        else:
            raise IllegalArgumentError('Invalid genomic region: {}'.format(region))

        # Update GR fields:
        self.region_str = region
        self.sites = self._region_str2sites()
        self.bp_tuple = (region_from, region_to)

    def _region_str2sites(self):
        # find CpG indexes in range of the region:
        # todo: find start and end separately (in case they are far apart)
        cmd = 'tabix {} {} | '.format(self.genome.dict_path, self.region_str)
        # cmd += 'awk \'{if (NR == 1) {first=substr($4,4)}}END{print first"-"substr($4,4)}\''
        cmd += 'awk \'{if (NR == 1) {first=$3}}END{print first"-"$3+1}\''
        # eprint(cmd)
        res = subprocess.check_output(cmd, shell=True).decode()

        # throw error if there are no CpGs in range
        if res.strip() == '-1':
            raise IllegalArgumentError('Invalid genomic region: {}. No CpGs in range'.format(self.region_str))

        s1, s2 = self._sites_str_to_tuple(res)
        # s2 += 1     # non-inclusive
        return s1, s2

    def _sites_str_to_tuple(self, sites_str):
        """ extract integers tuple (e.g (120, 130)) from a sites string (e.g '120-130') """
        if sites_str:
            sites_str = sites_str.replace(',', '')
            matchObj = re.match(r'([\d]+)-([\d]+)', sites_str)
            if matchObj:
                site1 = int(matchObj.group(1))
                site2 = int(matchObj.group(2))
                if not self.genome.nr_sites + 1 >= site2 > site1 >= 1:
                    msg = 'sites violate the constraints: '
                    msg += '{} >= {} > {} >= 1'.format(self.genome.nr_sites + 1, site2, site1)
                    raise IllegalArgumentError(msg)
                return site1, site2
        raise IllegalArgumentError('sites must be of format: ([\d])-([\d]).\nGot: {}'.format(sites_str))

    def index2chrom(self, site):
        if self.chrs_sz is None:
            self.chrs_sz = self.genome.get_chrom_cpg_size_table()
            self.chrs_sz['borders'] = np.cumsum(self.chrs_sz['size'])
        cind = np.searchsorted(np.array(self.chrs_sz['borders']).flatten(), site)
        return self.chrs_sz['chr'].loc[cind]

    def index2locus(self, index):
        """
        translate CpG index to genomic locus. e.g, CpG1 -> (chr1, 10469)
        :param index: a site index in range [1, NR_SITES]
        :return: chromosome, locus
        """
        index = int(index)
        # validate input
        if not self.genome.nr_sites + 1 >= index >= 1:
            eprint('Invalid site index:', index)
            raise IllegalArgumentError('Out of range site index:', index)

        # find locus:
        with open(self.genome.revdict_path, 'rb') as f:
            f.seek((index - 1) * 4)
            loc = np.fromfile(f, dtype=np.int32, count=1)[0] - 1
            return self.index2chrom(index), loc

    def __str__(self):
        if self.sites is None:
            return 'Whole genome'
        s1, s2 = self.sites
        f, t = self.bp_tuple
        # res = '{} ({:,} sites, {:,} bp, sites {}-{})'.format(self.region_str, s2 - s1, t - f + 1, s1, s2)
        res = '{} - {:,}bp, {:,}CpGs: {}-{}'.format(self.region_str, t - f + 1, s2 - s1, s1, s2)
        if self.annotation:
            res += '\n' + self.annotation
        return res

    def is_whole(self):
        """
        :return: True iff no filters (-r, -s) were applied. i.e, this gr is the whole genome.
        """
        return self.sites is None
