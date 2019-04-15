# view

View the content of input file (pat/unq/beta) as plain text.
Possible filter by genomic region or sites range
Output to stdout as default


Flags:
```
usage: view [-h] [-s SITES | -r REGION] [--genome GENOME] [-o OUT_PATH]
            [--sub_sample 0.0, 1.0] [-L BED_FILE] [--strict] [--inflate]
            [--awk_engine]
            input_file

View the content of input file (pat/unq/beta) as plain text. Possible filter
by genomic region or sites range Output to stdout as default

positional arguments:
  input_file

optional arguments:
  -h, --help            show this help message and exit
  -s SITES, --sites SITES
                        a CpG index range, of the form: "450000-450050"
  -r REGION, --region REGION
                        genomic region of the form "chr1:10,000-10,500"
  --genome GENOME       Genome reference name. Default is hg19.
  -o OUT_PATH, --out_path OUT_PATH
                        Output path. [stdout]
  --sub_sample (0.0, 1.0)
                        pat: subsample from reads.
  -L BED_FILE, --bed_file BED_FILE
                        pat: Only output reads overlapping the input BED FILE
  --strict              Truncate reads that start/end outside the given
                        region. Only relevant if "region", "sites" or
                        "bed_file" flags are given.
  --inflate             unq: add CpG-Index column to the output
  --awk_engine          pat: use awk engine instead of python. Its saves RAM
                        when dealing with large regions.
```

### random access (-r, -s flags)
View only reads (or values, in case of \*.beta files) overlapping the specified genomic region.
If no genomic region was provided, *view* outputs the whole file.
The genomic region may be specified in one of two ways:
1. CHROM:START-END, e.g. `-r chr1:10,747-10,758`
2. SITE1-SITE2, e.g. `-s 45-50`. This is non-inclusive, i.e. only sites 45,46,47,48,49 will be considered.
This feature is using *tabix* and the \*.csi index to achieve a quick random access (without reading the whole pat/unq file).

### inflate (unq)
the unq file does not store the CpG index as the pat file does.
Infer this field using the dictionary file "CpG.bed.gz", and output it along with the other columns.
```
% wgbs_tools view FILE.unq.gz -s 100000-100001
chr1	3649298	100	CCCCT	1
chr1	3649311	100	CCCTTC	1
% wgbs_tools view FILE.unq.gz -s 100000-100001 --inflate
chr1	99998	3649298	100	CCCCT	1
chr1	100000	3649311	100	CCCTTC	1
```

### awk_engine (pat)
view command loads the contents of the requested region into the RAM, which can be ineffective when specifying very large regions from large files. Using the *awk_engine* reduces memory usage and improves speed. However, it is a wrapper for awk, which might not be compatible with every systems.

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
