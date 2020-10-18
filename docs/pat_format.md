# pat.gz file

A gzipped tab separated text file with at least four columns: chromosome, CpG_index, methylation_pattern, and count. Other columns are optional.
Each line in the file stands for a read (molecule). 

* **chrom**: value is a string from _{chr1, chr2, …, chrX, chrY, chrM}_
* **CpG_index**: integer in range [1,NR_SITES]. The index of the first site occurring on the current read.
The pat file is sorted by this column. In [hg19](https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg19 "hg19 in UCSC"), NR_SITES == 28,217,448.


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


**Note:** The pat file is sorted by CpG_index column (`sort -k3,3n`), which is different from the UCSC order (`sort -k1,1 -k2,2n`)

#### compressing and indexing
The pat file is compressed using [bgzip](http://www.htslib.org/doc/bgzip.html) and indexed using [tabix](http://www.htslib.org/doc/tabix.html) (with a \*csi file). 
This compression is compatible with gzip, and with the [indexing](https://github.com/nloyfer/wgbs_tools/blob/master/docs/index.md), allows for a *fast random access* (see [view](https://github.com/nloyfer/wgbs_tools/blob/master/docs/view.md) command).


