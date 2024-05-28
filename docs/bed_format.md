# .bed files
If we want to use wgbstools with segments not found by the `segment` command, we need a specialized version of [`.bed`](https://en.wikipedia.org/wiki/BED_(file_format)) files.
The wgbstools `.bed` file has to have the following 5 columns: chr, start, end, startCpG, endCpG.  
Additional columns will not be used.  
We can use `convert` to add the required columns to a standard `.bed` file:

```bash
$ cat segments.bed
chr3	119527929	119528187
chr3	119528217	119528243
chr3	119528246	119528309
chr3	119528384	119528418
chr3	119528430	119528783

$ wgbstools convert -L segments.bed --out_path wgbs_segments.bed
$ cat wgbs_segments.bed
chr3	119527929	119528187	5394767	5394772	intron	NR1I2
chr3	119528217	119528243	5394774	5394777	intron	NR1I2
chr3	119528246	119528309	5394777	5394781	intron	NR1I2
chr3	119528384	119528418	5394782	5394786	intron	NR1I2
chr3	119528430	119528783	5394786	5394796	intron	NR1I2
```

As in other uses of convert, genomic annotations will be added when available (currently only hg19).

**Optional**: bgzip and index the '.bed' file, make it easier to access.
`index` wraps bgzip and tabix. It compresses a `.bed` (or `.pat`) file and generates a corresponding index file. **This step is necessary if you wish to visualize these blocks later** using the `vis` command.
```bash
$ wgbstools index wgbs_segments.bed
$ ls -1 wgbs_segments.*
wgbs_segments.bed.gz
wgbs_segments.bed.gz.tbi
```

- **TODO** - add reference to atlas segmentation
