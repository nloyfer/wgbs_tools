# Initialize a new genome reference


```
usage: init_genome [-h] [-f] genome_ref name

Init genome reference.

positional arguments:
  genome_ref   path to a genome *.fa file
  name         name of the genome (e.g. hg19, mm9...). A directory of this
               name will be created in references/

optional arguments:
  -h, --help   show this help message and exit
  -f, --force  Overwrite existing files if existed

```

For example,
```
wgbs_tools init_genome /path/to/my/genome.fa hg19
```
Will generate some necessary files in references/hg19/.
