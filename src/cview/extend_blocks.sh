#!/bin/bash
cut -f1,4,5 | sort -k2,2n | awk -v OFS='\t' '{print $1,(1>$2-100)?1:$2-100,$3}' | bedtools merge -i - -sorted
