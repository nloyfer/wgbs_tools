# beta file
beta file is the simplest representation of methylation data.
It is a binary file with a fixed size of 2 * 8 * NR_SITES bytes (~54MB for hg19),
holding a matrix of uint8 values with dimensions of (NR_SITES x 2).
For each of the NR_SITES CpG sites in the genome, it holds 2 values: the #meth and #covered. 
Meaning, the i'th row in the matrix corresponds to the i'th CpG site:
- **#meth**: the number of times the i'th site (CpGi) site is observed in a methylated state.
- **#coverage**: the total number of times i'th site (CpGi) is observed. #coverage==0 is equivalent to a missing value (NaN).

CpGi's beta value is obtained by dividing *#meth* / *#coverage*.

### reading a beta file, using
### python
```python
>>> import numpy as np
>>> content = np.fromfile(PATH, dtype=np.uint8).reshape((-1, 2))
>>> np.mean(content[:, 1])   # print average coverage
94.841405643770472
```
##### Fast Random Access with python:
```python
# reading only sites CpGi to CpGj (non-inclusive):
with open(PATH, 'rb') as f:
    f.seek((i - 1) * 2)
    res = np.fromfile(f, dtype=np.uint8, count=((j - i) * 2)).reshape((-1, 2))
```

### matlab
```matlab
fid=fopen(PATH,'r'); content=fread(fid); fclose(fid); A=reshape(content,2,[])';
```


### R
```R
> fname <- PATH
> N <- file.info(fname)$size
> content <- matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)
```

:exclamation:**Note**: this uint8 format limits values to range [0,255]. 
In case a CpG site appears over 255 times in the raw data (e.g. bam), its representation is normalized to this range. 
For example, if site CpG100 appeared 510 times, from which 100 times it was methylated, the 100'th row will be (50, 255).<br/>
We find it negligible for most cases, but if you find it important, use the `lbeta` format (`wgbstools pat2beta --lbeta ARGS...`)
