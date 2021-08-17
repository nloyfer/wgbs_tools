# view

View the content of input file (pat/beta) as plain text.
For beta files, view as bed.
Possible filter by genomic region or sites range
Output to stdout as default


Flags:
```
usage: view [-h] [--print_region] [-s SITES | -r REGION | -L BED_FILE]
            [--genome GENOME] [-o OUT_PATH] [--sub_sample [0.0, 1.0]]
            [--strict] [--strip] [-@ THREADS] [--min_len MIN_LEN]
            input_file

View the content of input file (pat/beta) as plain text. Possible filter by
genomic region or sites range Output to stdout as default

positional arguments:
  input_file

optional arguments:
  -h, --help            show this help message and exit
  --print_region        pat: Prints region before reads
  -s SITES, --sites SITES
                        a CpG index range, of the form: "450000-450050"
  -r REGION, --region REGION
                        genomic region of the form "chr1:10,000-10,500"
  -L BED_FILE, --bed_file BED_FILE
                        Bed file. Columns <chr, start, end>
  --genome GENOME       Genome reference name. Default is "default".
  -o OUT_PATH, --out_path OUT_PATH
                        Output path. [stdout]
  --sub_sample [0.0, 1.0]
                        pat: subsample from reads. Only supported for pat
  --strict              pat: Truncate reads that start/end outside the given
                        region. Only relevant if "region", "sites" or
                        "bed_file" flags are given.
  --strip               pat: Remove trailing dots (from beginning/end of
                        reads).
  -@ THREADS, --threads THREADS
                        Number of threads to use (default: all available CPUs)
  --min_len MIN_LEN     pat: Display only reads covering at least MIN_LEN CpG
                        sites [1]

```

### random access (-r, -s flags)
View only reads (or values, in case of \*.beta files) overlapping the specified genomic region.
If no genomic region was provided, *view* outputs the whole file.
The genomic region may be specified in one of two ways:
1. CHROM:START-END, e.g. `-r chr1:10,747-10,758`
2. SITE1-SITE2, e.g. `-s 45-50`. This is non-inclusive, i.e. only sites 45,46,47,48,49 will be considered.
This feature is using *tabix* and the \*.csi index to achieve a quick random access (without reading the whole pat file). For beta files, it utilizes the fact that they have fixed size, so random access takes O(1).

### strict (pat)
When specified with a region(s), *view* trims reads crossing the borders of the region. For example:
```
% wgbs_tools view FILE.pat.gz -s 167-168
chr1	166	TTTTT	1
chr1	167	TTTC	1
% wgbs_tools view FILE.pat.gz -s 167-168 --strict
chr1	167	T	1
chr1	167	T	1
```

### sub_sample (pat)
Subsamle from the pat file, *sub_sample* of the reads. the *count* field is taken into consideration.
For example, if a read has  *count*>1, it may be outputed with a smaller count.
```
% wgbs_tools view FILE.pat.gz -s 167-168
chr1	166	TTTTT	5
% wgbs_tools view FILE.pat.gz -s 167-168 --sub_sample .4
chr1	166	TTTTT	2
```
