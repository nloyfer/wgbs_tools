# wgbs_tools

#### Dependencies
- samtools
- python 3+


#### Quick start
First make sure you have a reference genome FASTA (e.g `hg19.fa`).
The FASTA must not be compressed (gzipped).
For exapmle, download it with:

```bash
curl https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -o hg19.fa.gz
gunzip hg19.fa.gz
```


```bash
# Download the repository
git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools

# compile the cpp files:
python3 setup.py

# Setup reference genome GENOME_NAME (e.g hg19).
python3 wgbs_tools.py init_genome /path/to/genome.fa GENOME_NAME
```

Now you can convert `pat.gz` and `beta` files out of `bam` files:
```bash
python3 wgbs_tools.py bam2pat BAM_PATH
```

#### Usage examples
Once you have `pat` and `beta` files, you can use wgbs_tools to visualize them. For example:

```bash
python3 wgbs_tools.py vis Sigmoid_Colon_STL003.pat.gz -r chr3:119528843-119529245
```
![alt text](https://github.com/nloyfer/wgbs_tools/blob/master/docs/img/pat_vis.png "pat vis example")

```bash
python3 wgbs_tools.py vis *.beta -r chr3:119528843-119529245
```
![alt text](https://github.com/nloyfer/wgbs_tools/blob/master/docs/img/beta_vis.png "beta vis example")
