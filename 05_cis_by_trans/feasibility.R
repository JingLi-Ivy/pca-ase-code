#!/usr/bin/env Rscript
# Feasibility: how many genes/samples are usable for the genome-wide ASE-based cis-trans screen?
suppressPackageStartupMessages({library(data.table)})
msg<-function(...)cat(sprintf(...),"\n")
setDTthreads(8)
a<-fread("/path/to/data_root/ASE/all_T_samples_SNPs_with_ASE_annotation.csv",
         select=c("gene_id","sample_id","ref.matches","alt.matches"))
a[,depth:=ref.matches+alt.matches]
a<-a[depth>0]
msg("raw records = %s ; genes = %s ; samples = %d", format(nrow(a),big.mark=","),
    format(uniqueN(a$gene_id),big.mark=","), uniqueN(a$sample_id))

gs<-a[,.(dep=sum(depth), ref=sum(ref.matches), nsnp=.N), by=.(gene_id,sample_id)]
msg("\n(gene,sample) pairs = %s", format(nrow(gs),big.mark=","))
cat("\n=== Distribution of total depth per (gene,sample) after aggregation ===\n"); print(summary(gs$dep))
print(quantile(gs$dep,c(.5,.75,.9,.95,.99)))
cat("\n=== Heterozygous sites per (gene,sample) ===\n"); print(summary(gs$nsnp))
msg("aggregation vs single top site: on average %.1f heterozygous sites per pair", mean(gs$nsnp))

for(D in c(8,10,20,30,50)){
  g2<-gs[dep>=D, .(nsamp=.N), by=gene_id]
  msg("\n--- threshold depth>=%d ---", D)
  msg("  eligible (gene,sample) pairs = %s", format(gs[dep>=D,.N],big.mark=","))
  for(S in c(20,30,50,80)) msg("  genes eligible in >=%d samples = %s", S, format(g2[nsamp>=S,.N],big.mark=","))
}
fwrite(gs,"/path/to/data_root/ASE/code/ase_gene_sample_agg.csv.gz")
msg("\nAggregated table saved: code/ase_gene_sample_agg.csv.gz")
