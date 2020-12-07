import re
import subprocess
import os
import os.path as op
import numpy as np
import pandas as pd
from io import StringIO
import multiprocessing
import sys

HG19_NR_SITES = 28217448  # total number of CpG sites in hg19 (fixed)
DIR = op.dirname(os.path.realpath(__file__)) + '/'

SRC_DIR = DIR + 'src/'
pat_sampler = SRC_DIR + 'pat_sampler/pat_sampler'
PAT2BETA_TOOL = SRC_DIR + 'pat2beta/stdin2beta'
PAT2RHO_TOOL = SRC_DIR + 'pat2beta/stdin2rho'
collapse_pat_script = SRC_DIR + 'collapse_pat.pl'
segment_tool = SRC_DIR + 'segment_betas/segmentor'
pat_segment_tool = SRC_DIR + 'segment_pats/pat_segmentor'

match_maker_tool = DIR + 'pipeline_wgbs/match_maker'
patter_tool = DIR + 'pipeline_wgbs/patter'

MAX_PAT_LEN = 150  # maximal read length in sites
MAX_READ_LEN = 1000  # maximal read length in bp

default_blocks_path = '/cs/cbio/netanel/blocks/outputs/blocks.s150.p15.bed.gz'  # todo: not generic

main_script = DIR + 'wgbs_tools.py'


class IllegalArgumentError(ValueError):
    pass


class GenomeRefPaths:
    def __init__(self, name='hg19'):
        self.genome = name
        self.refdir = self.build_dir()

        self.dict_path = self.join('CpG.bed.gz')
        self.chrom_cpg_sizes = self.join('CpG.chrome.size')
        self.chrom_sizes = self.join('chrome.size')
        self.revdict_path = self.join('CpG.rev.bin')
        self.genome_path = self.join('genome.fa')
        self.annotations = self.join('annotations.bed.gz', validate=False)
        if not op.isfile(self.annotations):
            self.annotations = None

        self.nr_sites = self.count_nr_sites()

    def count_nr_sites(self):
        if self.genome == 'hg19':
            return HG19_NR_SITES
        else:
            return int(self.get_chrom_cpg_size_table()['size'].sum())

    def join(self, fpath, validate=True):
        path = op.join(self.refdir, fpath)
        if validate and not op.isfile(path):
            raise IllegalArgumentError('Invalid reference path: ' + path)
        return path

    def build_dir(self):
        if not self.genome:
            self.genome = 'hg19'
        refdir = op.join(op.join(op.dirname(os.path.realpath(__file__)), 'references'), self.genome)
        if not op.isdir(refdir):
            raise IllegalArgumentError('Invalid reference name: {}'.format(self.genome))
        return refdir

    def get_chrom_cpg_size_table(self):
        return pd.read_csv(self.chrom_cpg_sizes, header=None, names=['chr', 'size'], sep='\t')


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class BedFileWrap:
    def __init__(self, bed_path, genome=None):
        self.bed_path = bed_path
        validate_single_file(bed_path)
        self.df = pd.read_csv(self.bed_path, usecols=[0, 1, 2], sep='\t',
                              names=['chr', 'start', 'end'], header=None, comment='#')

        # drop header line, if exists:
        if self.df.iloc[0, 0] == 'chr':
            self.df.drop(self.df.index[0], inplace=True)
            self.df.reset_index(inplace=True, drop=True)
        self.genome = genome

    def iter_grs(self):
        from genomic_region import GenomicRegion
        for _, r in self.df.iterrows():
            yield GenomicRegion(region='{}:{}-{}'.format(*r))
        # todo: check bed file has no overlaps?

    def fast_iter_regions(self):
        for _, r in self.df.iterrows():
            yield '{}:{}-{}'.format(*r)


def validate_dir(directory):
    if not op.isdir(directory):
        raise IllegalArgumentError('Invalid directory:\n{}'.format(directory))


def color_text(txt, cdict, scheme=16):
    if scheme == 16:
        return ''.join(['\033[{}m{}\033[00m'.format(cdict[c], c) if c in cdict.keys() else c for c in txt])
    elif scheme == 256:
        return ''.join(['\u001b[38;5;{}m{}\u001b[0m'.format(cdict[c], c) if c in cdict.keys() else c for c in txt])
    else:
        raise IllegalArgumentError('Invalid color scheme: {}'.format(scheme))


# def build_tool_path(tool_path):
#     return op.join(os.path.dirname(os.path.realpath(__file__)), tool_path)

def validate_prefix(prefix):
    """
    Return if prefix is a prefix of a file name inside a valid directory
    Raise exception otherwise
    """
    if op.isdir(prefix):
        raise IllegalArgumentError('Invalid prefix: {} is a directory.'.format(prefix))
    dirname = op.dirname(prefix)
    if not op.isdir(dirname):
        raise IllegalArgumentError('Invalid prefix: no such directory: {}'.format(dirname))


def load_dict(nrows=None, skiprows=None, genome_name='hg19'):
    d_path = GenomeRefPaths(genome_name).dict_path
    res = pd.read_csv(d_path, header=None, names=['chr', 'start'], sep='\t', usecols=[0, 1],
                      nrows=nrows, skiprows=skiprows)
    res['idx'] = res.index + 1
    return res


def load_dict_section(region, genome_name='hg19'):
    if region is None:
        return load_dict(genome_name = genome_name)
    cmd = 'tabix {} {}'.format(GenomeRefPaths(genome_name).dict_path, region)
    return read_shell(cmd, names=['chr', 'start', 'idx'])


def add_GR_args(parser, required=False, bed_file=False):
    """ Add genomic regions arguments parsing (flags [-s] and [-r]) """
    region_or_sites = parser.add_mutually_exclusive_group(required=required)
    region_or_sites.add_argument('-s', "--sites", help='a CpG index range, of the form: "450000-450050"')
    region_or_sites.add_argument('-r', "--region", help='genomic region of the form "chr1:10,000-10,500"')
    if bed_file:
        region_or_sites.add_argument('-L', "--bed_file", help='Bed file. Columns <chr, start, end>')
    parser.add_argument('--genome', help='Genome reference name. Default is hg19.', default='hg19')
    parser.add_argument('--no_anno', help='Do not print genomic annotations', action='store_true')


def add_multi_thread_args(parser):
    parser.add_argument('-@', '--threads', type=int, default=1,#multiprocessing.cpu_count(),
                        help='Number of threads to use (default: multiprocessing.cpu_count)')


def beta2vec(data, min_cov=1, na=np.nan):
    cond = data[:, 1] >= min_cov
    vec = np.divide(data[:, 0], data[:, 1], where=cond)  # normalize to range [0, 1)
    vec[~cond] = na
    return vec


def trim_to_uint8(data, lbeta = False):
    # Trim / normalize to range [0, 256)
    if lbeta:
        nr_bits = 16
        dtype = np.uint16
    else:
        nr_bits = 8
        dtype = np.uint8

    max_val = 2 ** nr_bits - 1
    big_indices = np.argwhere(data[:, 1] > max_val).flatten()
    data[:, 0][big_indices] = data[big_indices][:, 0] / data[big_indices][:, 1] * max_val
    data[:, 1][big_indices] = max_val
    return data.astype(dtype)


def load_dists(start, nr_sites, genome):
    """ load and return distance differences between adjacent sites """

    with open(genome.revdict_path, 'rb') as f:
        f.seek((start - 1) * 4)
        dists = np.fromfile(f, dtype=np.int32, count=nr_sites + 1)

    dists = dists[1:] - dists[:-1]
    dists[0] = 0
    dists[dists < 0] = 1e6 # cross chromosome hack
    return dists


def load_beta_data2(beta_path, gr=None, bed=None):
    if gr is not None and bed is not None:
        eprint('Error: both gr and bed_path supplied')
        raise IllegalArgumentError('Invalid usage of load_beta_data2')
    elif gr is not None and bed is None:
        return load_beta_data(beta_path, gr)
    elif gr is None and bed is not None:
        inds = load_dict_section(' -R ' + bed.bed_path, bed.genome)['idx'].values - 1
        return load_beta_data(beta_path)[inds, :]
    else:
        return load_beta_data(beta_path, None)



def load_beta_data(beta_path, sites=None):
    suff = op.splitext(beta_path)[1]
    if not (op.isfile(beta_path) and (suff in ('.beta', '.lbeta', '.bin'))):
        raise IllegalArgumentError("Invalid beta file:\n{}".format(beta_path))

    if suff == '.lbeta':
        sizet = 2
        dtype = np.uint16
    else:
        sizet = 1
        dtype = np.uint8

    if sites is None:
        data = np.fromfile(beta_path, dtype).reshape((-1, 2))
    else:
        start, end = sites
        with open(beta_path, 'rb') as f:
            f.seek((start - 1) * 2 * sizet)  # fix base-1 annotations
            data = np.fromfile(f, dtype=dtype, count=((end - start) * 2)).reshape((-1, 2))

    assert data.size, 'Data table is empty!'
    return data


def load_borders(borders_path, gr):
    validate_single_file(borders_path, '.gz')
    if not (op.isfile(borders_path + '.csi') or op.isfile(borders_path + '.tbi')):
        eprint('No csi found for file {}. Attempting to index it...'.format(borders_path))
        from index_wgbs import Indxer
        Indxer(borders_path).run()

    cmd = 'tabix {} {}'.format(borders_path, gr.region_str)
    res = subprocess.check_output(cmd, shell=True).decode()
    # print(cmd)
    df = pd.read_csv(StringIO(res), sep='\t', names=['start', 'end'], header=None, usecols=[3, 4])
    df = df - gr.sites[0]

    # remove blocks with < MIN_SITES_PER_BLOCKS sites:
    MIN_SITES_PER_BLOCKS = 1
    df = df[df['end'] - df['start'] >= MIN_SITES_PER_BLOCKS]

    # uniqe and sort all starts and ends:
    borders = set(df['start'])
    borders.update(set(df['end']))
    borders = np.array(sorted(borders))

    # keep only blocks within range [self.start, self.start + self.nsites]
    borders = borders[borders >= 0]
    borders = borders[borders <= gr.nr_sites]
    return borders


def validate_files_list(files, force_suff=None, min_len=1):
    """
    Make sure all files exist, and end with the same suffix.
    If force_suff is given, make sure it's identical (up to a leading '.') to the suffixes of the given files.
    If files are pat.gz or unq.gz, make sure their corresponding csi files exist.
    Make sure list has at least 2 items
    :param files: List of paths of files
    :param force_suff: None or an extension (e.g 'pat.gz', '.beta')
    :param min_len: integer. minimal number of items in the list
    :return: None. Raise an exception if the input is invalid.
    """

    if len(files) < min_len:
        raise IllegalArgumentError('Input error: at least {} input files must be given'.format(min_len))

    first = files[0]
    if len(first) == 1:
        raise IllegalArgumentError("Input is not a list of files:", files)

    # split extension of first file:
    suff = splitextgz(first)[1]

    if force_suff:
        if suff not in (force_suff, '.' + force_suff):
            raise IllegalArgumentError("Input files must be of type {}:".format(force_suff), files)

    # validate all files
    for file in files:
        validate_single_file(file, suff)


def validate_single_file(file, suff=None):
    """
    Make sure input file is valid:
        - file exists
        - has the correct suffix
        - has an index file, (for unq/pat). If not, attempt to create one.
    :param suff: 'beta', 'pat.gz' or 'unq.gz'
    """

    if file is None:
        raise IllegalArgumentError("Input file is None")

    if not (op.isfile(file)):
        raise IllegalArgumentError("No such file: {}".format(file))

    if not suff:
        suff = splitextgz(file)[1]
    elif not file.endswith(suff):
        raise IllegalArgumentError("file must end with {}:".format(suff), file)

    if suff in ('.pat.gz', '.unq.gz') and not op.isfile(file + '.csi'):
        eprint('No csi found for file {}. Attempting to index it...'.format(file))
        from index_wgbs import Indxer
        Indxer(file).run()


def splitextgz(input_file):
    """
    Extracts a file name + extension. If it's gz, add it to the other extension
    Examples:
        - fname.pat  -> (fname, pat)
        - fname.unq.gz -> (fname, unq.gz)
    """
    b, suff = op.splitext(input_file)
    if suff == '.gz':
        b, suff = op.splitext(b)
        suff += '.gz'
    return b, suff


def delete_or_skip(output_file, force):
    """
    If file exists, and force==True, remove it
    If file exists, and force==False, skip it, return False
    If file does not exist, return True.
    :param force: bool
    :return: False iff file should be skipped
    """
    # if file already exists, delete it or skip it
    if output_file is None:
        return True
    if op.isfile(output_file):
        if force:
            os.remove(output_file)
            if op.isfile(output_file + '.csi'):
                os.remove(output_file + '.csi')
        else:
            msg = 'File {} already exists. Skipping it.'.format(output_file)
            msg += 'Use [-f] flag to force overwrite.'
            eprint(msg)
            return False
    return True


def read_shell(command, **kwargs):
    """
    Based on timothymillar's code: https://github.com/pandas-dev/pandas/issues/16846
    Takes a shell command as a string and and reads the result into a Pandas DataFrame.

    Additional keyword arguments are passed through to pandas.read_csv.

    :param command: a shell command that returns tabular data
    :type command: str

    :return: a pandas dataframe
    :rtype: :class:`pandas.dataframe`
    """
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = proc.communicate()

    if proc.returncode == 0:
        txt = output.decode()
        if not txt:
            return pd.DataFrame()
        with StringIO(txt) as buffer:
            return pd.read_csv(buffer, sep='\t', header=None, **kwargs)
    else:
        message = ("Shell command returned non-zero exit status: {0}\n\n"
                   "Command was:\n{1}\n\n"
                   "Standard error was:\n{2}")
        raise IOError(message.format(proc.returncode, command, error.decode()))
