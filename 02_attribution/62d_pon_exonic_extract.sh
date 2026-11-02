#!/bin/bash
# 62d: S6 option A (chosen 2026-10-08): somatic variants after the panel-of-normals filter (filtered_with_pon_vcf, ANNOVAR-annotated;
# VAF >= 0.04, >= 5 alt reads, present in <= 1 PoN normal), ASE-cohort tumors only, no extra gnomAD filter (as Sen 2022).
# Streams each VCF and keeps exonic/splicing records: sample, chrom, pos, ref, alt, Func, Gene, ExonicFunc.
P=/path/to/data_root/variants/batch_45/mutect2_results/filtered_with_pon_vcf
O=/path/to/data_root/ASE/code/v8_round2/supp_data/S6_pon_exonic_variants.tsv
echo -e "sample\tchrom\tpos\tref\talt\tfunc\tgene\texonic_func" > $O
for t in /path/to/data_root/ASE/T*_ase_results.csv; do s=$(basename $t _ase_results.csv); id=${s#T}; f=$P/${id}_vs_${id}_filtered_pon.hg38_multianno.vcf
  [ -f $f ] || { echo "missing PoN VCF for $s" >&2; continue; }
  awk -F"\t" -v s=$s '!/^#/ && $8 ~ /Func.refGeneWithVer=(exonic|splicing|exonic\\x3bsplicing)/ {
     fu=$8; sub(/^(.*;)?Func.refGeneWithVer=/,"",fu); sub(/;.*/,"",fu); g=$8; sub(/.*Gene.refGeneWithVer=/,"",g); sub(/;.*/,"",g);
     ef=$8; sub(/.*ExonicFunc.refGeneWithVer=/,"",ef); sub(/;.*/,"",ef); print s"\t"$1"\t"$2"\t"$4"\t"$5"\t"fu"\t"g"\t"ef}' $f >> $O
done
echo "rows $(($(wc -l < $O)-1)) ; samples $(cut -f1 $O | sort -u | wc -l)"
