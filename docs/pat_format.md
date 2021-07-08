# pat.gz file

A bgzipped tab separated text file with at &ge;4 columns: chromosome, CpG_index, methylation_pattern, and count.<br/>
Additional columns are optional. Each line in the file stands for a read.

* **chrom**: value is a string from {`chr1`, `chr2`, ..., `chrX`, `chrY`, `chrM`}
* **CpG_index**: integer in range `[1,NR_SITES]`. The index of the first site occurring on the current read.
The pat file is sorted by this column. In [hg19](https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg19 "hg19 in UCSC"), `NR_SITES == 28,217,448`.

* **methylation pattern**: a string representing the methylation pattern of all consecutive CpG sites in the current read. 
Each site is represented by a single character: 
  * 'C' for methylated
  * 'T' for unmethylated
  * '.' for unknown state
* **count**: The number of times a read with these exact values (chrom, CpG_index and methylation pattern) occurred.

The rest of the nucleotides, and therefore the genomic distances between adjacent CpG sites, are omitted.

#### Example 
```
chr1	46	CCTC	2
chr1	47	CC...TC	1
chr1	47	T	13
```
The first read covers 4 CpG sites: CpG47 (methylated), CpG48 (methylated), CpG49 (*unmethylated*), and CpG50 (methylated). 
There are 2 reads with this exact pattern, starting at the same site.


:exclamation:**Note:** The pat file is sorted by CpG_index column (`sort -k2,2n -k3,3`. E.g., chr1 is followed by chr2, not chr10), which is different from the UCSC order (`sort -k1,1 -k2,2n`)

#### compressing and indexing
The pat file is compressed using [bgzip](http://www.htslib.org/doc/bgzip.html) and indexed using [tabix](http://www.htslib.org/doc/tabix.html) (with a \*csi file). 
This compression is compatible with gzip, and with the indexing allows for a *fast random access* (see [view](docs/view.md) command).


