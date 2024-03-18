
import argparse
import re
import shutil
import os
import os.path as op
import subprocess
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
from utils_wgbs import validate_single_file, eprint, IllegalArgumentError, \
        add_multi_thread_args, check_executable
from set_default_ref import set_def_ref


def generate_fai(fasta):
    """ Generate fai file if it does not exist """
    fai_path = fasta + '.fai'

    # If no fai file exists, return it
    if op.isfile(fai_path):
        return fai_path

    # Otherwise, generate it using samtools faidx
    eprint(f'[wt init] Indexing {fasta}')
    cmd = f'samtools faidx {fasta}'
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    output, error = p.communicate()
    # if failed to generate fai, print informative message and raise exception
    if p.returncode:
        eprint("[wt init] Failed with samtools idxstats %d\n%s\n%s" % (p.returncode, output.decode(), error.decode()))
        if fasta.endswith('.gz') and 'please use bgzip' in error.decode():
            msg = f'[wt init] Seems like your reference FASTA cannot be indexed with samtools faidx.\n' \
                    f'     Try one of the following:\n' \
                    f'     1. decompress it (gunzip {fasta}) and try again\n' \
                    f'     2. change the compression to bgzip:\n' \
                    f'        gunzip {fasta} && bgzip {fasta[:-3]}'
            eprint(msg)
        raise IllegalArgumentError('[wt init] Invalid reference FASTA')
    if op.isfile(fai_path):
        eprint(f'[wt init] Generated index file: {fai_path}')
    else:
        raise IllegalArgumentError('[wt init] Failed to generate index file (fai)')
    return fai_path


class InitGenome:
    def __init__(self, args):
        self.args = args
        self.ref_path = args.fasta_path
        self.force = args.force
        self.name = args.name

        # validate input files
        self.out_dir = self.setup_dir()
        self.get_fasta()
        self.fai_df = self.load_fai()

    def get_fasta(self):
        # download fasta from UCSC, unless the fasta file is provided
        if self.ref_path is not None:
            validate_single_file(self.ref_path)
            return

        # no FASTA path provided. Attempt to download one
        ref_path = op.join(self.out_dir, f'{self.name}.fa.gz')
        url = f'https://hgdownload.soe.ucsc.edu/goldenPath/{self.name}/bigZips/{self.name}.fa.gz'
        eprint(f'[wt init] No reference FASTA provided. Attempting to download from\n\t{url}')

        # make sure at least one of curl/ wget are installed
        if check_executable('curl'):
            cmd = f'curl {url} -o {ref_path}'
        elif check_executable('wget'):
            cmd = f'wget {url} -O {ref_path}'
        else:
            raise IllegalArgumentError('[wt init] Error: curl or wget must be installed')

        # Download
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        output, error = p.communicate()
        if p.returncode:
            eprint('[wt init] Failed downloading reference for genome ')
            eprint(f'{self.name}: %d\n%s' % (p.returncode, output.decode()))
            if error is not None:
                eprint(error.decode())
            eprint('[wt init] Try downloading yourself and use --fasta_name flag, or check the "name" parameter')
            raise IllegalArgumentError('[wt init] No reference FASTA found')
        eprint('[wt init] successfully downloaded FASTA. Now gunzip and bgzip it...')
        cmd = f'gunzip {ref_path} && bgzip -@ {self.args.threads} {ref_path[:-3]}'
        subprocess.check_call(cmd, shell=True)
        self.ref_path = ref_path

    def setup_dir(self):
        path = Path(op.realpath(__file__))
        out_dir = op.join(op.join(path.parent.parent.parent, 'references'),
                          self.name)
        # if the requested genome name already exist:
        # abort if --force was not specified, or delete the existing directory.
        if op.isdir(out_dir):
            if not self.force:
                msg = f'[wt init] Error: genome {self.name}'
                msg += f' already exists ({out_dir}).'
                msg += ' Use -f to overwrite it.'
                raise IllegalArgumentError(msg)
            else:
                shutil.rmtree(out_dir)
        os.makedirs(out_dir, exist_ok=True)
        eprint(f'[wt init] Setting up genome reference files in {out_dir}')
        return out_dir

    def load_fai(self):
        """ Generate, link and load the fai file to a DataFrame """
        fai_path = generate_fai(self.ref_path)

        # Link fa + fai (or fa.gz+fa.gz.gzi+fa.gz.fai) to the output dir
        fasta_name = 'genome.fa'
        if self.ref_path.endswith('.gz'):
            fasta_name += '.gz'
        self.link_file(self.ref_path, fasta_name)
        self.link_file(fai_path, fasta_name + '.fai')
        if fasta_name.endswith('.gz'):
            self.link_file(self.ref_path + '.gzi', fasta_name + '.gzi')

        # load fai file
        try:
            df = pd.read_csv(fai_path, sep='\t',
                             header=None,
                             usecols=[0, 1, 2, 4],
                             names=['chr', 'size', 'offset', 'width'])
            # filter invalid chromosomes
            df = df[df.apply(lambda x: is_valid_chrome(x['chr']), axis=1)]

            # sort chromosomes:
            if not self.args.no_sort:
                df = pd.DataFrame(sorted(df['chr'], key=chromosome_order),
                                  columns=['chr']).merge(df, how='left')
            return df
        except pd.errors.ParserError as e:
            raise IllegalArgumentError(f'Invalid fai file.\n{e}')

    def find_cpgs_loci(self):
        params = [(self.ref_path, c, self.args.debug)
                  for c in self.fai_df['chr']]
        p = Pool(self.args.threads)
        arr = p.starmap(load_seq_by_chrom, params)
        p.close()
        p.join()
        return pd.concat(arr).reset_index(drop=True)

    def run(self):
        eprint('[wt init] Processing chromosomes...')

        # CpG.bed.gz - Dictionary mapping locus to CpG-Index
        df = self.find_cpgs_loci()
        eprint('[wt init] Building CpG-Index dictionary...')
        df['site'] = df.index + 1
        # if not df.head()['chr'].unique()[0].startswith('chr'):
        #   df['chr'] = 'chr' + df['chr']
        self.dump_df(df, 'CpG.bed')
        cpg_dict = self.bgzip_tabix_dict(op.join(self.out_dir, 'CpG.bed'))

        # rev.CpG.bed.gz - symbolic link to the dictionary,
        # with an index file corresponds to the 3rd column, the CpG-Index,
        # To map CpG to locus
        rev_dict = op.join(self.out_dir, 'rev.CpG.bed.gz')
        os.symlink(cpg_dict, rev_dict)
        subprocess.check_call(f'tabix -b 3 -e 3 {rev_dict}', shell=True)

        # chrome.size - number of bases in each chromosome
        self.dump_df(self.fai_df[['chr', 'size']], 'chrome.size')

        # CpG.chrome.size - number of CpGs in each chromosome
        ncgs = self.fai_df[['chr']].copy()
        t = df.groupby('chr')['loc'].nunique().to_frame().\
            reset_index().rename(columns={'loc': 'size'})
        ncgs = ncgs.merge(t, on='chr', how='left').fillna(0)
        ncgs['size'] = ncgs['size'].astype(int)
        self.dump_df(ncgs[['chr', 'size']], 'CpG.chrome.size')

        # self.validate_nr_sites(df.shape[0])
        self.add_supp()  # add supplemental files for hg19
        eprint(f'[wt init] Finished initialization of genome {self.name}')

        # setup as default genome
        if not self.args.no_default:
            set_def_ref(self.name)

    def add_supp(self):
        path = Path(op.realpath(__file__))
        suppdir = op.join(path.parent.parent.parent, 'supplemental')

        if self.name == 'hg19':
            # link annotation files
            anno_file = op.join(suppdir, 'hg19.annotations.bed.gz')
            dst_anno = op.join(self.out_dir, 'annotations.bed.gz')
            if op.isfile(anno_file) and op.isfile(anno_file + '.tbi'):
                self.link_file(anno_file, dst_anno)
                self.link_file(anno_file + '.tbi', dst_anno + '.tbi')

            # link Illumina 450K file - hg19
            ilmn_file = op.join(suppdir, 'hg19.ilmn2CpG.tsv.gz')
            if op.isfile(ilmn_file):
                self.link_file(ilmn_file, op.join(self.out_dir, 'ilmn2CpG.tsv.gz'))

        elif self.name == 'hg38':
            # link Illumina 450K file
            ilmn_file = op.join(suppdir, 'hg38.ilmn2CpG.tsv.gz')
            if op.isfile(ilmn_file):
                self.link_file(ilmn_file, op.join(self.out_dir, 'ilmn2CpG.tsv.gz'))

    def validate_nr_sites(self, nr_sites):
        if self.args.debug:
            return
        d = {
            'mm9': 13120864,
            'hg19': 28217448
        }
        if self.name in d.keys():
            if nr_sites != d[self.name]:
                msg = '[wt init] WARNING: number of sites of the reference genome '
                msg += f'{self.name} is usually {d[self.name]}, but you got {nr_sites}'
                eprint(msg)

    def link_file(self, src, dst):
        if not op.isfile(src):
            raise IllegalArgumentError(f'[wt init] Invalid reference genome file: {src}')
        src = op.abspath(src)
        dst = op.join(self.out_dir, dst)
        cmd = f'ln -s {src} {dst}'
        if op.islink(dst):
            os.unlink(dst)
        subprocess.check_call(cmd, shell=True)

    def bgzip_tabix_dict(self, dict_path):
        eprint('[wt init] bgzip and index...')
        subprocess.check_call(f'bgzip -@ {self.args.threads} -f {dict_path}', shell=True)
        subprocess.check_call(f'tabix -Cf -b 2 -e 2 {dict_path}.gz', shell=True)
        return dict_path + '.gz'

    def dump_df(self, df, path):
        path = op.join(self.out_dir, path)
        df.to_csv(path, index=None, header=False, sep='\t')


def load_seq_by_chrom(ref_path, chrom, debug=False):
    eprint(f'[wt init] chromosome: {chrom}')

    # load the chromosome's subsequence from fasta
    cmd = f'samtools faidx {ref_path} {chrom} | tail -n +2'
    if debug:
        cmd += ' | head -200'
    txt = subprocess.check_output(cmd, shell=True).decode()
    seq = ''.join(s.strip() for s in txt.split('\n')).upper()

    # Find CpG sites loci
    tf = pd.DataFrame([m.start() + 1 for m in re.finditer('CG', seq)],
                      columns=['loc'])
    tf['chr'] = chrom
    return tf[['chr', 'loc']]


def chromosome_order(c):
    if c.startswith('chr'):
        # raise IllegalArgumentError('Invalid chromosome' + c)
        c = c[3:]
    if c.isdigit():
        return int(c)
    if c == 'X':
        return 10000
    if c == 'Y':
        return 10001
    if c in ('M', 'MT'):
        return 10002
    return 10003


def is_valid_chrome(chrome):
    # A chromosome is valid if it has the form "chrX",
    # where X is digit(s) or (X,Y,M)
    return bool(re.match(r'^(chr)?([\d]+|[XYM]|(MT))$', chrome))


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('name', help='name of the genome (e.g. hg19, mm9...).\n'
                                     'A directory of this name will be created '
                                     'in references/')
    parser.add_argument('--fasta_path',
            help='path to a reference genome FASTA file.\n' \
                 ' If none provided, wgbstools will attempt to download one from UCSC')
    parser.add_argument('-f', '--force', action='store_true',
            help='Overwrite existing files if existed')
    parser.add_argument('--no_default', action='store_true',
            help='Do not set this genome as the default reference for wgbstools.\n' \
                    'This setting can be changed with wgbstools set_default_ref')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--no_sort', action='store_true',
                        help='If set, keep the chromosome order of the reference genome.\n'
                             'Default behaviour is to sort 1,2,...,10,11,...,X,Y,M')
    add_multi_thread_args(parser)
    args = parser.parse_args()
    return args


def main():
    """ Init genome reference. """
    args = parse_args()
    InitGenome(args).run()


if __name__ == '__main__':
    main()
