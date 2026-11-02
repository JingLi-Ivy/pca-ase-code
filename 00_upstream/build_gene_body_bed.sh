#!/usr/bin/env bash
# Gene-body BED (chrom, start0, end, Ensembl gene_id; chromosome names without "chr") from the
# Ensembl GRCh38 release 110 GTF. Read by run_star_wasp_pipeline_{tumor,normal}.sh to map SNPs to genes.
# Requires gawk (three-argument match()).
set -euo pipefail

GTF=/path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz
OUT=/path/to/data_root/wasp/gtf_gene_body_geneid.bed

zcat "$GTF" \
| awk 'BEGIN{OFS="\t"} $3=="gene"{chr=$1; gsub(/^chr/,"",chr); start=$4-1; if(start<0) start=0; end=$5; m=match($0,/gene_id "([^"]+)"/,a); gid=(m?a[1]:"NA"); print chr,start,end,gid}' \
| sort -k1,1V -k2,2n \
> "$OUT"

echo "DONE: $OUT"
