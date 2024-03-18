import re
import subprocess
import os
import os.path as op
from pathlib import Path
import numpy as np
from utils_wgbs import IllegalArgumentError, GenomeRefPaths, eprint


def index2chrom(site, genome):
    chrs_sz = genome.get_chrom_cpg_size_table()
    return chrs_sz['chr'].iloc[np.searchsorted(np.cumsum(chrs_sz['size'].values), site)]


def get_genome_name(gname):
    if gname is None or gname == 'default':
        path = Path(op.realpath(__file__))
        refdir = op.join(op.join(path.parent.parent.parent, 'references'), 'default')
        return os.readlink(refdir)
    return gname


class GenomicRegion:
    def __init__(self, args=None, region=None, sites=None, array_id=None, genome_name=None):
        self.genome_name = get_genome_name(genome_name)
        self.chrom = None
        self.sites = sites
        self.array_id = array_id
        self.region_str = region
        self.bp_tuple = None
        self.args = args

        # todo: this could be prettier
        if args is not None:
            self.genome_name = get_genome_name(args.genome)
            self.genome = GenomeRefPaths(self.genome_name)
            if args.sites:
                self.parse_sites(args.sites)
            elif args.region:
                self.parse_region(args.region)
            elif args.array_id:
                self.parse_array_id(args.array_id)
        elif region is not None:
            self.genome = GenomeRefPaths(self.genome_name)
            self.parse_region(region)
        elif sites is not None:
            self.genome = GenomeRefPaths(self.genome_name)
            self.parse_sites(sites)
        elif array_id is not None:
            self.genome = GenomeRefPaths(self.genome_name)
            self.parse_array_id(array_id)
        else:
            raise IllegalArgumentError(f'Invalid GR init {region}')

        self.nr_sites = None if self.sites is None else self.sites[1] - self.sites[0]
        self.annotation = self.add_anno()

    def add_anno(self):
        if self.args is None or self.is_whole() or 'no_anno' not in self.args:
            return
        if self.args.no_anno:
            return
        anno_path = self.genome.annotations
        if anno_path is None:
            return
        try:
            cmd = f'tabix {anno_path} {self.region_str} | cut -f4- | uniq'
            return subprocess.check_output(cmd, shell=True).decode().strip()
        except subprocess.CalledProcessError:
            eprint(f'Failed to retrieve annotation for reagion {self.region_str}')

    def parse_sites(self, sites_str):
        """ Parse input of the type -s / --sites (e.g 15-25) """

        # Parse sites string:
        s1, s2 = self._sites_str_to_tuple(sites_str)

        # Translate sites indexes to genomic loci:
        self.chrom, region_from = self.index2locus(s1)
        chrom2, region_to = self.index2locus(s2 - 1)  # non-inclusive
        region_to += 1  # include the whole last site (C and G)
        if self.chrom != chrom2:
            eprint(f'ERROR: sites range cross chromosomes! ({s1}, {s2})')
            raise IllegalArgumentError('Invalid sites input')

        # Update GR fields:
        self.sites = (s1, s2)
        self.region_str = f'{self.chrom}:{region_from}-{region_to}'
        self.bp_tuple = (region_from, region_to)

    def _chrome_size(self):
        df = self.genome.get_chrom_size_table()
        return int(df[df['chr'] == self.chrom]['size'].values[0])

    def find_region_format(self, region):
        region = region.replace(',', '')  # remove commas

        # In case region is a whole chromosome
        chrome_match = re.match(r'^(chr)?([\d]+|[XYM]|(MT))$', region)
        if chrome_match:
            if region not in self.genome.get_chroms():
                raise IllegalArgumentError(f'Unknown chromosome: {region}')
            self.chrom = region
            return region, 1, self._chrome_size()

        # match region string to format chrom:from
        uni_region_match = re.match(r'^(chr)?([\d]+|[XYM]|(MT)):([\d]+)$', region)
        if uni_region_match:
            region_from = uni_region_match.group(4)
            region += f'-{int(region_from) + 1}'

        # match region string to format chrom:from-to
        region_match = re.match(r'^((chr)?([\d]+|[XYM]|(MT))):([\d]+)-([\d]+)$', region)
        if not region_match:
            raise IllegalArgumentError(f'Invalid genomic region: {region}')

        self.chrom = region_match.group(1)
        if self.chrom not in self.genome.get_chroms():
            raise IllegalArgumentError(f'Unknown chromosome: {region}')
        region_from = int(region_match.group(5))
        region_to = int(region_match.group(6))

        return region, region_from, region_to


    def parse_region(self, region):
        """ Parse input of the type -r / --region (e.g chr11:200-300) """

        self.region_str, region_from, region_to = self.find_region_format(region)

        # validate region range:
        if region_to <= region_from:
            raise IllegalArgumentError(f'Invalid genomic region: {region}. end before start')
        if region_to > self._chrome_size() or region_from < 1:
            raise IllegalArgumentError(f'Invalid genomic region: {region}. Out of range')

        # Update GR fields:
        self.bp_tuple = (region_from, region_to)
        self.sites = self._region_str2sites()

    def _region_str2sites(self):
        # find CpG indexes in range of the region:
        cmd = f'tabix {self.genome.dict_path} {self.region_str} | '
        # cmd += 'awk \'(NR==1){first=$3} {lbp=$2} END{print first"-"$3+1}\''

        # if bp_tuple[1] equals exactly a loci of a CpG site, this site is *not* included
        # e.g., in hg19, chr6:71046415-71046562 is 9718430-9718435 in sites
        cmd += f"awk -v b={self.bp_tuple[1]} "
        cmd += '\'(NR==1){first=$3} END{if ($2<b) {r+=1}; print first"-"$3+r}\''
        res = subprocess.check_output(cmd, shell=True).decode()
        # eprint(cmd)

        if len(set(res.strip().split('-'))) == 1:
            res = '-1'

        # throw error if there are no CpGs in range
        if res.strip() == '-1':
            raise IllegalArgumentError(f'Invalid genomic region: {self.region_str}. No CpGs in range')

        s1, s2 = self._sites_str_to_tuple(res)
        return s1, s2

    def _sites_str_to_tuple(self, sites_str):
        """ extract integers tuple (e.g (120, 130)) from a sites string (e.g '120-130') """
        if not sites_str:
            raise IllegalArgumentError(f'Empty sites string: {sites_str}')

        sites_str = sites_str.replace(',', '')
        # start-end syntax
        matchObj = re.match(r'([\d]+)-([\d]+)', sites_str)
        if matchObj:
            site1 = int(matchObj.group(1))
            site2 = int(matchObj.group(2))
        # single site syntax:
        elif '-' not in sites_str and sites_str.isdigit():
            site1 = int(sites_str)
            site2 = site1 + 1
        else:
            raise IllegalArgumentError(f'sites must be of format: "start-end" or "site" .\nGot: {sites_str}')
        # validate sites are in range:
        if not self.genome.get_nr_sites() + 1 >= site2 >= site1 >= 1:
            msg = 'sites violate the constraints: '
            msg += f'{self.genome.get_nr_sites() + 1} >= {site2} > {site1} >= 1'
            raise IllegalArgumentError(msg)
        if site1 == site2:
            site2 += 1
        return site1, site2

    def index2locus(self, index):
        """
        translate CpG index to genomic locus. e.g, CpG1 -> (chr1, 10469)
        :param index: a site index in range [1, NR_SITES]
        :return: chromosome, locus
        """
        index = int(index)
        # validate input
        if not self.genome.get_nr_sites() + 1 >= index >= 1:
            eprint('Invalid site index:', index)
            raise IllegalArgumentError('Out of range site index:', index)

        # find chromosome:
        chrom = index2chrom(index, self.genome)
        # find locus:
        cmd = f'tabix {self.genome.revdict_path} {chrom}:{index}-{index} | cut -f2'
        try:
            loc = int(subprocess.check_output(cmd, shell=True).decode().strip())
        except ValueError as e:
            msg = f'Failed retrieving locus for site {index} with command:\n{cmd}\n{e}'
            raise IllegalArgumentError(msg)
        return chrom, loc

    def parse_array_id(self, array_id):
        """ Parse input of the type --array_id (e.g cg00001755) """

        # validate ID
        if not (array_id.startswith('cg') and len(array_id) > 2 and array_id[2:].isdigit()):
            eprint(f'ERROR: Invalid Illumina array id: {array_id}')
            raise IllegalArgumentError('Invalid Illumina array ID')
        # verify there is an Illumina map file
        idict = self.genome.ilmn2cpg_dict
        if idict is None or not op.isfile(idict):
            raise IllegalArgumentError(f'Could not find Illumina map file: {idict}')

        # Find cg ID in the map file
        try:
            cmd = f'gunzip -c {idict} | grep -w {array_id} | cut -f2'
            cpg_ind = int(subprocess.check_output(cmd, shell=True).decode().strip())
        except ValueError as e:
            msg = f'Failed retrieving locus for site {array_id} with command:\n{cmd}\n{e}'
            raise IllegalArgumentError(msg)

        self.parse_sites(str(cpg_ind))

    def __str__(self):
        if self.sites is None:
            return 'Whole genome'
        s1, s2 = self.sites
        nr_bp = np.diff(self.bp_tuple)[0] + 1
        res = f'{self.region_str} - {nr_bp:,}bp, {s2 - s1:,}CpGs: {s1}-{s2}'
        if self.annotation:
            res += '\n' + self.annotation
        return res

    def is_whole(self):
        """ True iff no filters (-r, -s) were applied.
            i.e, this gr is the whole genome."""
        return self.sites is None
