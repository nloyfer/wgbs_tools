import subprocess
import os
import os.path as op
from io import StringIO
import multiprocessing
import sys
from pathlib import Path
import shutil
import numpy as np
import pandas as pd


current_file_path = Path(op.realpath(__file__))
DIR = str(current_file_path.parent)

SRC_DIR = op.join(current_file_path.parent.parent.parent, 'src/')
pat_sampler = SRC_DIR + 'pat_sampler/pat_sampler'
pat2beta_tool = SRC_DIR + 'pat2beta/stdin2beta'
mask_pat_tool = SRC_DIR + 'pat2beta/mask_pat'
collapse_pat_script = SRC_DIR + 'collapse_pat.pl'
segment_tool = SRC_DIR + 'segment_betas/segmentor'
cview_tool = SRC_DIR + 'cview/cview'
cview_extend_blocks_script = SRC_DIR + 'cview/extend_blocks.sh'
view_beta_script = SRC_DIR + 'view_beta.sh'
view_lbeta_script = SRC_DIR + 'view_lbeta.sh'
homog_tool = SRC_DIR + 'homog/homog'
add_loci_tool = SRC_DIR + 'cpg2bed/add_loci'
plot_marker_heatmap_script = SRC_DIR + 'R/plot_marker_heatmap.R'
plot_tree_script = SRC_DIR + 'R/plot_circle.R'

match_maker_tool = SRC_DIR + 'pipeline_wgbs/match_maker'
patter_tool = SRC_DIR + 'pipeline_wgbs/patter'
add_cpg_count_tool = SRC_DIR + 'pipeline_wgbs/add_cpg_counts'
allele_split_tool = SRC_DIR + 'pipeline_wgbs/snp_patter'
bam_meth_split_tool = SRC_DIR + 'pipeline_wgbs/bam_split.sh'


MAX_PAT_LEN = 150  # maximal read length in sites

COORDS_COLS3 = ['chr', 'start', 'end']
COORDS_COLS5 = COORDS_COLS3 + ['startCpG', 'endCpG']


main_script = op.join(DIR, 'wgbs_tools.py')


class IllegalArgumentError(ValueError):
    pass

class EmptyBamError(IllegalArgumentError):
    pass

class GenomeRefPaths:
    def __init__(self, name=None):
        self.genome = name
        self.refdir = self.build_dir()

        self.dict_path = self.join('CpG.bed.gz')
        self.chrom_cpg_sizes = self.join('CpG.chrome.size')
        self.chrom_sizes = self.join('chrome.size')
        self.revdict_path = self.join('rev.CpG.bed.gz')
        self.genome_path = self.join('genome.fa', validate=False)
        self.annotations = self.join('annotations.bed.gz', validate=False)
        self.blocks = self.join('blocks.bed.gz', validate=False)
        self.blacklist = self.join('blacklist.bed', validate=False)
        self.whitelist = self.join('whitelist.bed', validate=False)
        self.ilmn2cpg_dict = self.join('ilmn2CpG.tsv.gz', validate=False)
        self._chrom_cpg_size_table = None
        self._chrome_size_table = None
        self._chroms = None

    def get_nr_sites(self):
        return int(self.get_chrom_cpg_size_table()['size'].sum())

    def join(self, fpath, validate=True):
        path = op.join(self.refdir, fpath)
        if not op.isfile(path):
            if op.isfile(path + '.gz'):
                path += '.gz'
            else:
                if validate:
                    raise IllegalArgumentError('Invalid reference path: ' + path)
                path = None
        return path

    def build_dir(self):
        if not self.genome:
            self.genome = 'default'
        path = Path(op.realpath(__file__))
        refdir = op.join(op.join(path.parent.parent.parent, 'references'), self.genome)
        if self.genome == 'default':
            self.genome = os.readlink(refdir)
            refdir = str(Path(refdir).resolve())

        if not op.isdir(refdir):
            raise IllegalArgumentError(f'Invalid reference name: {self.genome}')
        return refdir

    def get_chrom_cpg_size_table(self):
        if self._chrom_cpg_size_table is None:
            self._chrom_cpg_size_table = pd.read_csv(self.chrom_cpg_sizes,
                    header=None, names=['chr', 'size'], sep='\t')
        return self._chrom_cpg_size_table

    def get_chrom_size_table(self):
        if self._chrome_size_table is None:
            self._chrome_size_table = pd.read_csv(self.chrom_sizes,
                    header=None, names=['chr', 'size'], sep='\t')
        return self._chrome_size_table

    def get_chroms(self):
        if self._chroms is None:
            self._chroms = tuple(pd.read_csv(self.chrom_sizes,
                    header=None, names=['chr'], usecols=[0], sep='\t')['chr'].values)
        return self._chroms


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def mkdirp(dpath):
    # mkdir "dpath" if not exist or None
    if dpath:
        Path(dpath).mkdir(parents=True, exist_ok=True)
    return dpath


def check_samtools_version(major=1, minor=15, verbose=False):
    # make sure samtools version >= major.minor
    try:
        t = subprocess.check_output(['samtools', '--version'])
        existing_version = t.splitlines()[0].strip().split()[1].decode()
        if verbose:
            eprint('[wt] samtools path:', shutil.which('samtools'), sep='\t')
            eprint('[wt] samtools version:', existing_version, sep='\t')
        nums = existing_version.split('.')
        cmajor = int(nums[0])
        cminor = int(nums[1])
        if cmajor != major:
            return cmajor > major
        return cminor >= minor
    except Exception:
        eprint('[wt] WARNING: failed to validate samtools version')
        # failed to run samtools --version
        pass
    return False


def check_executable(cmd, verbose=False):
    for p in os.environ['PATH'].split(":"):
        if os.access(op.join(p, cmd), os.X_OK):
            return True
    if verbose:
        eprint(f'executable {cmd} not found in PATH')
    return False


def validate_local_exe(tool):
    # make sure the file exists
    if not op.isfile(tool):
        eprint(f'[wt] binary executable not found: {tool}. Run setup.py to compile')
        raise IllegalArgumentError('Missing executable')
    # make sure it's executable
    if not os.access(tool, os.X_OK):
        eprint(f'[wt] file {tool} is not executable')
        raise IllegalArgumentError('Invalid executable')


def validate_out_dir(out_dir, verbose=True):
    if not out_dir:
        out_dir = '.'
    if not op.isdir(out_dir):
        if verbose:
            eprint(f'[wt] creating output directory {out_dir}')
        os.mkdir(out_dir)
    if not os.access(out_dir, os.W_OK | os.X_OK):
        raise IllegalArgumentError('Output directory has no writing permissions')


def drop_dup_keep_order(lst):
    seen = set()
    return [x for x in lst if not (x in seen or seen.add(x))]


def validate_dir(directory):
    if not op.isdir(directory):
        raise IllegalArgumentError(f'Invalid directory:\n{directory}')
    return directory


def color_text(txt, cdict, scheme=16):
    if scheme not in [16, 256]:
        raise IllegalArgumentError(f'Invalid color scheme: {scheme}')
    if scheme == 16:
        res = ''.join([f'\033[{cdict[c]}m{c}\033[00m' if c in cdict.keys() else c for c in txt])
    elif scheme == 256:
        res = ''.join([f'\u001b[38;5;{cdict[c]}m{c}\u001b[0m' if c in cdict.keys() else c for c in txt])

    return res

# def build_tool_path(tool_path):
#     return op.join(os.path.dirname(os.path.realpath(__file__)), tool_path)

def validate_prefix(prefix):
    """
    Return if prefix is a prefix of a file name inside a valid directory
    Raise exception otherwise
    """
    if op.isdir(prefix):
        raise IllegalArgumentError(f'Invalid prefix: {prefix} is a directory.')
    dirname = op.dirname(prefix)
    if not op.isdir(dirname):
        raise IllegalArgumentError(f'Invalid prefix: no such directory: {dirname}')


def load_dict(nrows=None, skiprows=None, genome_name=None):
    d_path = GenomeRefPaths(genome_name).dict_path
    res = pd.read_csv(d_path, header=None, names=['chr', 'start'], sep='\t', usecols=[0, 1],
                      nrows=nrows, skiprows=skiprows)
    res['idx'] = res.index + 1
    return res


def load_dict_section(region, genome_name=None):
    if region is None:
        return load_dict(genome_name=genome_name)
    dpath = GenomeRefPaths(genome_name).dict_path
    cmd = f'tabix {dpath} {region}'
    return read_shell(cmd, names=['chr', 'start', 'idx'])


def add_GR_args(parser, required=False, bed_file=False, no_anno=False, expand=False):
    """ Add genomic regions arguments parsing (flags [-s] and [-r]) """
    region_or_sites = parser.add_mutually_exclusive_group(required=required)
    region_or_sites.add_argument('-s', '--sites', help='a CpG index range, of the form: "450000-450050"')
    region_or_sites.add_argument('-r', '--region', help='genomic region of the form "chr1:10,000-10,500"')
    region_or_sites.add_argument('--array_id', help='Illumina array id, e.g. cg00001755')
    if bed_file:
        region_or_sites.add_argument('-L', '--bed_file', help='Bed file. Columns <chr, start, end>. '\
                'For some features columns 4-5 should be <startCpG, endCpG> (run wgbstools convert -L BED_PATH)')
    parser.add_argument('--genome', help='Genome reference name. Default is "default".', default='default')
    if no_anno:
        parser.add_argument('--no_anno', help='Do not print genomic annotations', action='store_true')
    if expand: # todo: implement, use in vis.
        parser.add_argument('--expand', help='expand region by K bp or by K sites for each side', type=int)
    return region_or_sites


def add_multi_thread_args(parser):
    try:
        cpu_env = 'SLURM_JOB_CPUS_PER_NODE'
        if cpu_env in os.environ.keys():
            def_cpus = int(os.environ[cpu_env])
        else:
            def_cpus = multiprocessing.cpu_count()
    except Exception:
        def_cpus = 8
    parser.add_argument('-@', '--threads', type=int, default=def_cpus,
                        help='Number of threads to use (default: all available CPUs)')


def add_no_beta_arg(parser):
    parser.add_argument('--no_beta', action='store_true', help='Do not generate a beta file')

def add_no_pat_arg(parser):
    parser.add_argument('--no_pat', action='store_true', help='Do not generate a pat file')


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


def beta_sanity_check(beta_path, genome):
    # sanity test: make sure beta file has the correct number of sites
    # (fits current genome)
    nr_sites_in_beta = op.getsize(beta_path) // 2
    if beta_path.endswith('.lbeta'):
        nr_sites_in_beta /= 2
    if int(nr_sites_in_beta) != genome.get_nr_sites():
        eprint(f'[wt beta] WARNING: beta file size ({nr_sites_in_beta:,} sites)\n' \
               f'          incomatible with current genome reference ' \
               f'({genome.get_nr_sites():,} sites)')
        return False
    return True


def load_beta_data(beta_path, sites=None, genome=None):
    if genome is not None:
        beta_sanity_check(beta_path, genome)
    suff = op.splitext(beta_path)[1]
    if not (op.isfile(beta_path) and (suff in ('.beta', '.lbeta', '.bin'))):
        raise IllegalArgumentError(f'Invalid beta file:\n{beta_path}')

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

    assert data.size, beta_path + ': Data table is empty!'
    return data


def load_borders(bpath, gr, genome):
    if bpath is False:
        return np.array([])
    if bpath is True:
        bpath = GenomeRefPaths(genome).blocks
        if bpath is None:
            eprint(f'[wt blocks] default blocks path not found: {bpath}')
            return np.array([])

    # else, bpath is a string

    validate_single_file(bpath, '.bed.gz')
    if not op.isfile(bpath + '.tbi'):
        eprint(f'No tbi found for file {bpath}. Attempting to index it...')
        from index import Indxer
        Indxer(bpath).run()

    df = read_shell(f'tabix {bpath} {gr.region_str}', usecols=[3, 4])     # load borders section
    borders = np.sort(np.unique(df.values.flatten())) - gr.sites[0]       # sort, unique, shift
    return borders[np.logical_and(borders >= 0, gr.nr_sites >= borders)]  # return only blocks in range


def validate_file_list(files, force_suff=None, min_len=1):
    """
    Make sure all files exist, and end with the same suffix.
    If force_suff is given (e.g. "beta"), make sure all files ends with it
    If files are pat.gz make sure their corresponding csi files exist.
    Make sure list has at least min_len items
    :param files: List of paths of files
    :param force_suff: None or an extension (e.g 'pat.gz', '.beta')
    :param min_len: integer. minimal number of items in the list
    :return: None. Raise an exception if the input is invalid.
    """

    if len(files) < min_len:
        raise IllegalArgumentError(f'Input error: at least {min_len} input files must be given')

    first = files[0]
    if len(first) == 1:
        raise IllegalArgumentError(f'Input is not a list of files: {files}')

    if (force_suff is not None) and (not first.endswith(force_suff)):
        raise IllegalArgumentError(f'Input file {first} must end with {force_suff}')

    # validate all files
    suff = splitextgz(first)[1]
    for fpath in files:
        validate_single_file(fpath, suff)


def validate_single_file(fpath, suff=None):
    """
    Make sure input fpath is valid:
        - fpath exists
        - has the correct suffix
        - has an index fpath, (for pat). If not, attempt to create one.
    :param suff: 'beta' or 'pat.gz'
    """

    if fpath is None:
        raise IllegalArgumentError("Input file is None")

    if not op.isfile(fpath):
        raise IllegalArgumentError(f'No such file: {fpath}')

    if suff is not None and not fpath.endswith(suff):
        raise IllegalArgumentError(f'file {fpath} must end with {suff}')

    if fpath.endswith('.pat.gz') and not op.isfile(fpath + '.csi'):
        eprint(f'No csi found for file {fpath}. Attempting to index it...')
        from index import Indxer
        Indxer(fpath).run()

    return fpath


def splitextgz(input_file):
    """
    Extracts a file name + extension. If it's gz, add it to the other extension
    Examples:
        - fname.pat  -> (fname, pat)
        - fname.pat.gz -> (fname, pat.gz)
    """
    b, suff = op.splitext(input_file)
    if suff == '.gz':
        b, suff = op.splitext(b)
        suff += '.gz'
    return b, suff


def pretty_name(fpath):
    return splitextgz(op.basename(fpath))[0]


def safe_remove(fpath):
    if fpath is not None and op.isfile(fpath):
        os.remove(fpath)

def mult_safe_remove(files):
    for f in files:
        safe_remove(f)

def delete_or_skip(output_file, force):
    """
    If file exists, and force==True, remove it
    If file exists, and force==False, skip it, return False
    If file does not exist, return True.
    :param force: bool
    :return: False iff file should be skipped
    """
    # if file already exists, delete it or skip it
    if output_file is None or output_file == sys.stdout or output_file == '/dev/stdout':
        return True
    if op.isfile(output_file):
        if force:
            mult_safe_remove([output_file, output_file + '.csi'])
        else:
            msg = f'File {output_file} already exists. Skipping it.'
            msg += ' Use [-f] flag to force overwrite.'
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


def bed2reg(df):
    if not set(COORDS_COLS3).issubset(set(df.columns)):
        raise IllegalArgumentError('[wt] missing coordinate columns in bed file')
    return df['chr'].astype(str) + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)


def catch_BrokenPipeError():
    os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
    sys.exit(1)
