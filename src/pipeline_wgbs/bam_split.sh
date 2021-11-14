#!/bin/bash
homog_prop=$1;
min_cpg=$2;
awk -v homog_prop="$homog_prop" -v min_cpg="$min_cpg" '{
	for (i=1;i<=NF;i++) {
		if ($i ~ /YI/) {
			split($i,methyl_info,":");
			split(methyl_info[3],methyl_counts,",");
			total_cpgs = methyl_counts[1] + methyl_counts[2];
			if (total_cpgs >= min_cpg){
				unmethyl_prop = methyl_counts[2] / total_cpgs;
				methyl_prop = methyl_counts[1] / total_cpgs;
                if (homog_prop >= 0.5){
                    if (methyl_prop >= homog_prop){
                        print $0;
                    }
                } else {
                    if (methyl_prop <= homog_prop){
                        print $0;
                    }
                }
			}
	    }
	}
}'
