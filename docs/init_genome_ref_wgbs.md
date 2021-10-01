# Initialize a new genome reference

Setup a reference genome (e.g., hg19, mm9). This is a necessary step before using any other command.
This command generates a directory with the genome name in the "references" directory, and saves there about 100-200Mb of indexing files that are necessary for other wgbstools commands.
When you setup a reference genome, it becomes the default genome for wgbstools from now on, until you set another genome as the default (`wgbstools set_default_ref`).

```
usage: init_genome [-h] [--fasta_path FASTA_PATH] [-f] [--no_default] [-d]
                   [--no_sort] [-@ THREADS]
                   name

Init genome reference. Note: we currently support only chromosomes starting
with "chr". I.e. "chr1" and not "1".

positional arguments:
  name                  name of the genome (e.g. hg19, mm9...). A directory of
                        this name will be created in references/

optional arguments:
  -h, --help            show this help message and exit
  --fasta_path FASTA_PATH
                        path to a reference genome FASTA file. If none
                        provided, wgbstools will attempt to download one from
                        UCSC
  -f, --force           Overwrite existing files if existed
  --no_default          Do not set this genome as the default reference for
                        wgbstools. This setting can be changed with wgbstools
                        set_default_ref
  -d, --debug
  --no_sort             If set, keep the chromosome order of the reference
                        genome. Default behaviour is to sort
                        1,2,...,10,11,...,X,Y,M
  -@ THREADS, --threads THREADS
                        Number of threads to use (default: all available CPUs)
```

For example,
```
wgbs_tools init_genome hg19 --fasta_path /path/to/my/genome.fa
```
Will generate some necessary files in `references/hg19/`.

If no `fasta_path` is specified, `wgbstools` will attempt to download a fasta file with that name from UCSC.
The default reference genome for any wgbstools command from now on will be `hg19`.

#### Multiple reference genomes
wgbstools supports multiple reference genomes. The default genome is the last genome configured using the `init_genome` command. To change the default genome, run `set_default_ref`, for example:
```bash
# download and setup the mm9 reference with:
wgbs_tools init_genome mm9
# Now the default ref is mm9. change it back to hg19 with:
wgbstools set_default_ref --name hg19
# to list the configured genomes and the default one, run:
wgbstools set_default_ref -ls
```
