#!/bin/bash
dict=$1
beta=$2
outpath=$3
paste <(gunzip -c $dict | awk -v OFS='\t' '{print $1,$2,$2+2}') <(hexdump -v -e '/1 "%u\t" /1 "%u\n"' $beta) > $outpath
