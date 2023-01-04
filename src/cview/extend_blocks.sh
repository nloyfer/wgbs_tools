#!/bin/bash

cut -f1,4,5 | sort -k2,2n | grep -v -w "NA" | grep -v '^$' | \
awk -v OFS="\t" '
{
    start = $2 - 100;
    if (start < 1) {
        start = 1;
    }
    print $1,start,$3
}
' | bedtools merge -i - -sorted
