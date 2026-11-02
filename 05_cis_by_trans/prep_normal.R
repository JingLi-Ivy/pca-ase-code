#!/usr/bin/env Rscript
# Prepare the T/N-specificity test: normal-side aggregated depth + a composition proxy comparable across tissues (marker score) + proxy validation
suppressPackageStartupMessages({library(data.table)}); setDTthreads(16)
msg<-function(...){cat(sprintf(...),"\n"); flush.console()}
OUT<-"/path/to/data_root/ASE/code/cistrans_full"; dir.create(OUT,showWarnings=FALSE)

cat("########## 1. Building the SNP->gene map from the tumor annotation ##########\n")
ann<-fread("/path/to/data_root/ASE/all_T_samples_SNPs_with_ASE_annotation.csv",
           select=c("seqtnames","start","gene_id"))
map<-unique(ann); rm(ann); gc()
map[,seqtnames:=as.character(seqtnames)]
msg("unique (site,gene) mappings = %s ; unique sites = %s",
    format(nrow(map),big.mark=","), format(uniqueN(map[,.(seqtnames,start)]),big.mark=","))
setkey(map,seqtnames,start)

cat("\n########## 2. Aggregated depth in normal samples ##########\n")
fs<-list.files("/path/to/data_root/ASE",pattern="^N\\d+_allele_specific_counts\\.pipeline\\.nonzero\\.txt$",full.names=TRUE)
msg("normal sample files = %d", length(fs))
agN<-rbindlist(lapply(fs,function(f){
  x<-fread(f,header=FALSE,select=c(1,2,6,7),col.names=c("seqtnames","start","ref","alt"))
  x[,seqtnames:=as.character(seqtnames)]
  y<-map[x,on=.(seqtnames,start),nomatch=0]
  y[,.(dep=sum(ref+alt), nsnp=.N, ref=sum(ref)), by=gene_id][,sid:=gsub("\\D","",basename(f))][]
}),fill=TRUE)
msg("normal (gene,sample) pairs = %s", format(nrow(agN),big.mark=","))
msg("normal aggregated depth: median=%.0f  mean heterozygous sites=%.1f", median(agN$dep), mean(agN$nsnp))
fwrite(agN,file.path(OUT,"ase_gene_sample_agg_N.csv.gz"))

cat("\n########## 3. Composition proxy comparable across tissues (marker score) ##########\n")
fpkm<-fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls")
mk<-list(
 epi=c("EPCAM","KRT8","KRT18","KRT5","KRT15","CDH1","KLK3","NKX3-1"),
 immune=c("PTPRC","CD3E","CD2","CD19","MS4A1","CD68","LYZ","ITGAM","CD53","LCP1"),
 stroma=c("COL1A1","COL1A2","DCN","LUM","ACTA2","PDGFRB","THY1"))
score<-function(cols,genes){
  M<-as.matrix(fpkm[gene_name %in% genes,..cols]); if(nrow(M)==0) return(rep(NA,length(cols)))
  Z<-t(scale(t(log2(M+0.01)))); colMeans(Z,na.rm=TRUE)}
res<-list()
for(tis in c("T","N")){
  cc<-grep(paste0("^",tis,"\\d+_WTS$"),names(fpkm),value=TRUE)
  res[[tis]]<-data.table(sid=gsub("\\D","",cc), tissue=tis,
    epi_s=score(cc,mk$epi), imm_s=score(cc,mk$immune), str_s=score(cc,mk$stroma))}
sc<-rbindlist(res)
msg("Composition score: T %d, N %d", nrow(res$T), nrow(res$N))
fwrite(sc,file.path(OUT,"composition_scores_TN.csv"))

cat("\n########## 4. Proxy validation: against the true deconvolution fractions on the tumor side ##########\n")
cf<-fread("/path/to/data_root/singlecell/scRNA_deconvolution/GSE278936/deconvolution_T/output/CPGEA_cell_fraction.tsv")
nm<-names(cf)[1]; cf[,sid:=gsub("\\D","",get(nm))]
num<-setdiff(names(cf)[sapply(cf,is.numeric)],"sid")
ep<-grep("pitheli",num,value=TRUE)[1]
imm<-grep("T_|B_|NK|Monocyte|Macrophage|DC|Mast|Neutrophil",num,value=TRUE)
cf[,`:=`(epi_d=rowSums(.SD)),.SDcols=ep][,imm_d:=rowSums(.SD),.SDcols=imm]
v<-merge(res$T, cf[,.(sid,epi_d,imm_d)], by="sid")
msg("comparable tumor samples = %d", nrow(v))
msg("Spearman(marker epithelial score, deconvolved epithelial fraction) = %+.3f", cor(v$epi_s,v$epi_d,method="spearman"))
msg("Spearman(marker immune score, deconvolved immune fraction) = %+.3f", cor(v$imm_s,v$imm_d,method="spearman"))
msg("Spearman(marker stromal score, deconvolved epithelial fraction) = %+.3f (expected negative)", cor(v$str_s,v$epi_d,method="spearman"))
msg("-> if both epithelial and immune rho are >0.7, the proxy can be used equivalently in both tissues")
cat("\nComposition score difference between tissues (T vs N):\n")
for(k in c("epi_s","imm_s","str_s"))
  msg("  %s : T median=%+.3f  N median=%+.3f  Wilcoxon p=%.3g", k,
      median(res$T[[k]],na.rm=TRUE), median(res$N[[k]],na.rm=TRUE),
      wilcox.test(res$T[[k]],res$N[[k]])$p.value)
