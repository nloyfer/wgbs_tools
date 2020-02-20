# Initialize a new genome reference


```
usage: init_genome [-h] [-f] [-d] [--no_sort] [-@ THREADS] genome_ref name

Init genome reference.

positional arguments:
  genome_ref            path to a genome *.fa file
  name                  name of the genome (e.g. hg19, mm9...). A directory of
                        this name will be created in references/

optional arguments:
  -h, --help            show this help message and exit
  -f, --force           Overwrite existing files if existed
  -d, --debug
  --no_sort             If set, keep the chromosome order of the reference
                        genome. Default behaviour is to sort
                        1,2,...,10,11,...,X,Y,M
  -@ THREADS, --threads THREADS
                        Number of threads to use (default:
                        multiprocessing.cpu_count)

```

For example,
```
wgbs_tools init_genome /path/to/my/genome.fa hg19
```
Will generate some necessary files in references/hg19/.
