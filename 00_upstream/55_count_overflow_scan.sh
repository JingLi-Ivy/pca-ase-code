#!/bin/bash
# 55: scan all per-sample WASP allele count files for 16-bit overflow (negative ref/alt counts, cols 6-7; values near 32767).
# Output: one line per affected site x sample, and per-sample summary.
cd /path/to/data_root/ASE
OUT=/path/to/data_root/ASE/code/v8_round2/count_overflow; mkdir -p $OUT
: > $OUT/overflow_sites.tsv
for f in [TN]*_allele_specific_counts.pipeline.nonzero.txt; do s=${f%%_*}
  awk -v s=$s '($6<0 || $7<0 || $6>=30000 || $7>=30000){print s"\t"$1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7}' $f >> $OUT/overflow_sites.tsv; done
awk -F"\t" '{neg=($6<0||$7<0); n[$1]++; m[$1]+=neg} END{for(k in n) print k"\t"n[k]"\t"m[k]}' $OUT/overflow_sites.tsv | sort > $OUT/per_sample.tsv
echo "rows $(wc -l < $OUT/overflow_sites.tsv); negative rows $(awk -F"\t" '$6<0||$7<0' $OUT/overflow_sites.tsv | wc -l); samples $(wc -l < $OUT/per_sample.tsv)"
