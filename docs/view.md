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
                        pat: subsample from reads. Only supported for pat
  -L BED_FILE, --bed_file BED_FILE
                        pat: Only output reads overlapping the input BED FILE
  --strict              Truncate reads that start/end outside the given
                        region. Only relevant if "region", "sites" or
                        "bed_file" flags are given.
  --inflate             unq: add CpG-Index column to the output
  --awk_engine          pat: use awk engine instead of python. Its saves RAM
                        when dealing with large regions.
```

### inflate
the unq file does not store the CpG index as the pat file does.
Infers this field using the dictionary file "CpG.bed.gz", and output it.
```
% wgbs_tools view FILE.unq.gz -s 100000-100001
chr1	3649298	100	CCCCT	1
chr1	3649311	100	CCCTTC	1
% wgbs_tools view FILE.unq.gz -s 100000-100001 --inflate
chr1	99998	3649298	100	CCCCT	1
chr1	100000	3649311	100	CCCTTC	1
```
