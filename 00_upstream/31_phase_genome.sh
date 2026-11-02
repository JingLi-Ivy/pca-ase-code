#!/bin/bash
# Genome-wide statistical phasing of CPGEA germline SNPs (Beagle 5.5 + 1000G GRCh38 reference), for Flynn-style aFC.
# Read-only inputs (combined.134.qc, 1000G); outputs under v8_round2/phase_genome/. Each step is skipped if its output exists;
# the 1000G subset is deleted only after phasing succeeded.
# absolute paths: xargs subshells do not inherit the conda PATH
BIN=$HOME/.conda/envs/r_env/bin
O=/path/to/data_root/ASE/code/v8_round2/phase_genome; mkdir -p $O; cd $O
PP=/path/to/data_root/CPGEA/CPGEA_genotype/combined.134.qc
G=/path/to/data_root/reference/variation/1000G/1000G_GRCh38_2019
BJ=/path/to/data_root/ASE/code/v8_round2/tools/beagle.27Feb25.75f.jar
one() {
  c=$1; [ -s phased_chr$c.vcf.gz ] && { echo "chr$c already phased"; return; }
  R=$G/ALL.chr$c.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
  [ -s cpgea_chr$c.vcf.gz ] || $BIN/plink2 --bfile $PP --chr $c --snps-only just-acgt --max-alleles 2 --export vcf bgz --out cpgea_chr$c --memory 8000 --threads 2 > /dev/null
  [ -s ref_chr$c.vcf.gz.csi ] || { $BIN/bcftools view -v snps -m2 -M2 $R -Oz -o ref_chr$c.vcf.gz && $BIN/bcftools index -f ref_chr$c.vcf.gz; }
  [ -s refal_chr$c.tsv ] || $BIN/bcftools query -f "%POS\t%REF\t%ALT\n" ref_chr$c.vcf.gz > refal_chr$c.tsv
  [ -s aligned_chr$c.vcf.gz.csi ] || { $BIN/Rscript /path/to/data_root/ASE/code/v8_round2/31a_align_fast.R cpgea_chr$c.vcf.gz refal_chr$c.tsv aligned_chr$c.vcf > align_chr$c.log 2>&1 \
    && $BIN/bcftools view aligned_chr$c.vcf -Oz -o aligned_chr$c.vcf.gz && $BIN/bcftools index -f aligned_chr$c.vcf.gz && rm -f aligned_chr$c.vcf; }
  [ -s aligned_chr$c.vcf.gz.csi ] || { echo "chr$c FAILED at alignment"; return; }
  java -Xmx16g -jar $BJ gt=aligned_chr$c.vcf.gz ref=ref_chr$c.vcf.gz out=phased_chr$c impute=false nthreads=4 seed=1 > beagle_chr$c.out 2>&1
  if [ -s phased_chr$c.vcf.gz ]; then $BIN/bcftools index -f phased_chr$c.vcf.gz; rm -f ref_chr$c.vcf.gz ref_chr$c.vcf.gz.csi refal_chr$c.tsv; echo "chr$c phased $(date)"
  else echo "chr$c FAILED at Beagle (see beagle_chr$c.out)"; fi
}
export -f one; export PP G BJ BIN
seq 1 22 | xargs -P 4 -I{} bash -c 'one {}'
echo "phased chromosomes: $(ls phased_chr*.vcf.gz 2>/dev/null | wc -l) / 22"
