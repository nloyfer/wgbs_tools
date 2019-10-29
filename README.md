# wgbs_tools

## Quick start
```bash
# Download the repository
git clone URL
cd wgbs_tools
# Setup reference genome names GENOME_NAME (e.g hg19).
# This step may take a few minutes.
python3 wgbs_tools.py init_genome /path/to/genome.fa GENOME_NAME

compile the cpp files:
g++ -std=c++11 src/pat2beta/stdin2beta.cpp -o src/pat2beta/stdin2beta
g++ -std=c++11 src/pat_sampler/sampler.cpp -o src/pat_sampler/pat_sample
g++ -std=c++11 wgbs_tools/pipeline_wgbs/patter.cpp -o wgbs_tools/pipeline_wgbs/patter
g++ -std=c++11 wgbs_tools/pipeline_wgbs/match_maker.cpp -o wgbs_tools/pipeline_wgbs/match_maker

```
