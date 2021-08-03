# wgbstools tutorial
## installation and configuration
First install `wgbstools` and configure the `hg19` genome
```bash
git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools
python setup.py
wgbstools init_genome hg19
```

It is recommended to add wgbstools to your $PATH, E.g,
```bash
export PATH=${PATH}:$PWD
```

## all set! Let's begin
### Data and region
For this short tutorial, we will use the following publicly available samples from ROADMAP:

| SRX  | Tissue  |  Donor |
|---|---|---|
| SRX175350 |  Lung cells          | STL002
| SRX388743 |  Pancreas cells      | STL002
| SRX190161 |  Sigmoid colon cells | STL003

```bash
$ cd tutorial
$ ls -1 *bam
Lung_STL002.small.bam
Pancreas_STL002.small.bam
Sigmoid_Colon_STL003.small.bam
```

To keep things compact, we consider a small region of ~4Kb, covering 100 CpG sites.
`convert` command translates genomic loci to CpG-index range and vice verca. It also prints genomic annotations, when available (currently only hg19).
```bash
$ region=chr3:119527929-119531943
$ wgbstools convert -r $region
chr3:119527929-119531943 - 4,015bp, 100CpGs: 5394767-5394867
intron  NR1I2
exon    NR1I2
intron  NR1I2
exon    NR1I2
```
### Generate pat & beta files
To generage `pat` and `beta` files for each of the 3 samples, we use the `bam2pat` command.
```bash
$ wgbstools bam2pat *.bam -r $region
```

### Segmentation 
Segment the region into homogenously methylated blocks
```bash
$ wgbstools segment --betas *beta --min_cpg 3 --max_bp 2000 -o blocks.small.bed -r $region
```
bgzip and index the output blocks, make it easier to access
```bash
$ wgbstools index blocks.small.bed
```

### collapse beta files to blocks
Collapse the beta files to the blocks we just found:
```bash
$ wgbstools beta_to_table blocks.small.bed.gz --betas *beta
```

### Visualizations
try different forms of visualizations
```bash
$ wgbstools vis -r chr3:119527929-119531943 -b blocks.small.bed.gz *beta
$ wgbstools vis -r chr3:119527929-119531943 -b blocks.small.bed.gz *beta --heatmap
$ wgbstools vis -r chr3:119528585-119528783 -b blocks.small.bed.gz Sigmoid_Colon_STL003.small.pat.gz --min_len 4
$ wgbstools vis -r chr3:119528585-119528783 -b blocks.small.bed.gz Sigmoid_Colon_STL003.small.pat.gz --min_len 4 --strict
$ wgbstools view -s 5394796-5394834 Sigmoid_Colon_STL003.small.pat.gz --sub_sample .05
```
