# unq.gz file

The unq format is an extension for the [pat format](https://github.com/nloyfer/wgbs_tools/blob/master/docs/pat_format.md).
This format also keeps the information about reads lengths, and higher resolution locations.
ItA gzipped tab separated file with five columns: chromosome, start, length, methylation_pattern, and count.
chrom, methylation_pattern and count are similar as in the pat format.
* **start**: The read's beginning locus.
* **length**: The read's length in base-pair.
The CpG_index field is missing, since it's easy to infer it from the *start* field.

#### Example 
```
chr1	10712	100	CCTC	1
chr1	10713	100	CCTC	1
chr1	10724	100	CC...TC	1
chr1	10751	100	T	2
chr1	10752	100	T	2
```

