# Index

bgzip and index pat or unq files. Accepts single or multiple files. 
Files may be in various states: bgzipped, gzipped or uncompressed (i.e extensions .pat[.gz] or .unq[.gz]).
bgzips them and generate an index for them (csi).


```
usage: index [-h] [-f] input_files [input_files ...]

positional arguments:
  input_files  One or more file with extensions .pat[.gz] or .unq[.gz]

optional arguments:
  -h, --help   show this help message and exit
  -f, --force  Overwrite existing index file (csi) if existed

```
