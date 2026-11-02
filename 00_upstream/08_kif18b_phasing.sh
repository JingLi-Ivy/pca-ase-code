#!/bin/bash
# A7 + A2: statistical phasing of rs9914422 and the KIF18B ASE reporter SNPs (Beagle 5.5 + 1000G GRCh38 reference) -> phase-resolved, DNA-corrected allelic ratios
exec > /path/to/data_root/ASE/code/v8_round2/08_kif18b_phasing.log 2>&1
set -eo pipefail
source /path/to/software/share/conda3/etc/profile.d/conda.sh; conda activate r_env
O=/path/to/data_root/ASE/code/v8_round2/kif18b_phase; mkdir -p $O; cd $O
A=/path/to/data_root/ASE/code/v8_round2/kif18b_audit
REF=/path/to/data_root/reference/variation/1000G/1000G_GRCh38_2019/ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools view -r 17:44870000-45000000 -v snps -m2 -M2 $REF -Oz -o ref1kg.vcf.gz; bcftools index -f ref1kg.vcf.gz
bcftools query -f "%POS\t%REF\t%ALT\n" ref1kg.vcf.gz > ref_alleles.tsv
Rscript /path/to/data_root/ASE/code/v8_round2/08a_align_alleles.R $A/kif18b_region.vcf.gz ref_alleles.tsv cpgea_aligned.vcf
bcftools view cpgea_aligned.vcf -Oz -o cpgea_aligned.vcf.gz; bcftools index -f cpgea_aligned.vcf.gz
java -Xmx8g -jar /path/to/data_root/ASE/code/v8_round2/tools/beagle.27Feb25.75f.jar gt=cpgea_aligned.vcf.gz ref=ref1kg.vcf.gz out=phased impute=false nthreads=4 seed=1
bcftools index -f phased.vcf.gz
Rscript /path/to/data_root/ASE/code/v8_round2/08b_phased_allelic_ratio.R
echo "exit: $?"
