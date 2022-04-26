#!/bin/bash
dict=$1
lbeta=$2
if [ "$#" -eq 2 ]; then
    paste <(gunzip -c $dict | awk -v OFS='\t' '{print $1,$2-1,$2+1}') <(hexdump -v -e '/2 "%u\t" /2 "%u\n"' $lbeta)
elif [ "$#" -eq 5 ]; then
    chrom=$3
    start=$4
    nr_sites=$5
    end=$(( start+nr_sites-1 ))
    paste <(tabix $dict ${chrom}":"${start}"-"${end} | awk -v OFS='\t' '{print $1,$2-1,$2+1}') <(hexdump -s $(( 4*(start-1) )) -n $(( 4*nr_sites )) -v -e '/2 "%u\t" /2 "%u\n"' $lbeta)
else
    echo "Illegal number of parameters"
fi

