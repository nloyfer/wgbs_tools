# pat.gz file

A bgzipped tab-separated text file storing the read-level methylation information.
#### Example:
```
chr1    46       CT    1
chr1    47       CC..TC  1
chr1    47       T       13
chr2    2300000  C       4
chr10   14633440 TTC     1
```
The first read covers 2 CpG sites: CpG46 (methylated) and CpG47 (unmethylated). 

#### Specification
A pat files contains &ge;4 columns:<br/>

* **chrom**: a string (E.g., `chr1`, `chr2`, ..., `chrX`, `chrY`, `chrM`).
* **CpG index**: integer in range `[1,NR_SITES]`. The index of the first site occurring on the current read.
The pat file is sorted by this column. <br/>
In [hg19](https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg19 "hg19 in UCSC"), `NR_SITES == 28,217,448`. Run `wgbstools convert -s CpG_i-CpG_j` to translate the CpG indexes to genomic loci.

* **methylation pattern**: a string representing the methylation pattern of all consecutive CpG sites in the current read. 
Each site is represented by a single character: <br/>

| symbol  | meaning  |
|---|---|
| `C`  | Methylated CpG  |
| `T`  | Unmethylated CpG  |
| `.`  | CpG with an unknown status  |


* **count**: The number of times a read with these exact values (chrom, CpG index and methylation pattern) occurred.
* Additional columns are optional. 

### From bam to pat
A pat file can be thought of as a compact methylation-oriented representation of a bam file. 
Most of the fields are omitted, as well as non-CpG nucleotides. Reads not covering any CpG sites are discarded.
In paired-end data, the two mates are merged into a single read in the pat file.
When 2+ reads cover the exact same CpG sites and present identical methylation pattern, they are merged to a single line in the pat file, and the count field is increased. The rest of the nucleotides are omitted.


:exclamation:**Note:** The pat file is sorted by CpG index column (`sort -k2,2n -k3,3`. E.g., chr1 is followed by chr2, not by chr10), which is different from the UCSC order (`sort -k1,1 -k2,2n`)

#### compressing and indexing
The pat file is compressed using [bgzip](http://www.htslib.org/doc/bgzip.html) and indexed using [tabix](http://www.htslib.org/doc/tabix.html) (with a `*.csi` file). 
This compression is compatible with gzip, and with the indexing allows for a *fast random access* (try `wgbstools view`).


