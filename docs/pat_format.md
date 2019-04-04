# pat.gz file

A gzipped tab separated file with four columns: chromosome, CpG_index, methylation_pattern, and count.
Each line in the file stands for a read (molecule). 

* **chrom**: value is a string from _{chr1, chr2, …, chrX, chrY, chrM}_
* **CpG_index**: integer in range [1,NR_SITES]. The index of the first site occurring on the current read.
The pat file is sorted by this column. In [hg19](https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg19 "hg19 in UCSC"), NR_SITES == 28217448.



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
There are 2 reads with exact pattern, starting at the same site.


:exclamation:**Note:** The pat file is sorted by CpG_index column, which is different from the genome browser order (sort -k1,1 -k2,2n)

#### remarks / future work:


