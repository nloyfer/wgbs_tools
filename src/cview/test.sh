#!/usr/bin/zsh
# test w/o flags
echo "test no flags"
cat <(cut -f4,5 t.bed) <(echo -1) <(zcat p.pat.gz ) | cview | sort -k2,2n -k3,3 > c.tsv
wgbstools view -L t.bed  p.pat.gz | sort -k2,2n -k3,3 -u > p.tsv
diff p.tsv c.tsv | wcl
rm p.tsv c.tsv

# test --strict
echo "--strict"
cat <(cut -f4,5 t.bed) <(echo -1) <(zcat p.pat.gz ) | cview --strict | sort -k2,2n -k3,3 > c_strict.tsv
wgbstools view -L t.bed --strict p.pat.gz | sort -k2,2n -k3,3 > p_strict.tsv
diff c_strict.tsv p_strict.tsv | wcl
rm c_strict.tsv p_strict.tsv

# test --strict --strip
echo "--strict --strip"
cat <(cut -f4,5 t.bed) <(echo -1) <(zcat p.pat.gz ) | cview --strict --strip | sort -k2,2n -k3,3 > c_strict_strip.tsv
wgbstools view -L t.bed --strict --strip p.pat.gz | sort -k2,2n -k3,3 > p_strict_strip.tsv
diff c_strict_strip.tsv p_strict_strip.tsv | wcl
rm c_strict_strip.tsv p_strict_strip.tsv

# test --strict --strip --min
echo "--strict --strip --min_len"
cat <(cut -f4,5 t.bed) <(echo -1) <(zcat p.pat.gz ) | cview --strict --strip --min_cpgs 3 | sort -k2,2n -k3,3 > c_strict_strip_min_cpg.tsv
wgbstools view -L t.bed --strict --strip --min_len 3 p.pat.gz | sort -k2,2n -k3,3 > p_strict_strip_min_cpg.tsv
diff p_strict_strip_min_cpg.tsv c_strict_strip_min_cpg.tsv
rm p_strict_strip_min_cpg.tsv c_strict_strip_min_cpg.tsv
