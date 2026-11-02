#!/usr/bin/env Rscript
# TRUE copy-neutral test: per-sample mean allelic imbalance restricted to BALANCED DIPLOID
# regions (FACETS tcn=2, lcn=1) -> purely regulatory imbalance. Does it still rise with Gleason?
suppressPackageStartupMessages({ library(data.table) })
PROJ <- "/path/to/data_root/ASE"; OUT <- file.path(PROJ,"figures/ase_clinical")
FACETS <- "/path/to/data_root/variants/batch_43/facets/results"
SPLIT  <- file.path(PROJ,"allelic_direction/ase_split")
num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]","",as.character(x))))
lg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), sprintf(...)))

fits <- list.files(FACETS, pattern="_facets_fit\\.rds$", full.names=TRUE)
samples <- sub("_facets_fit\\.rds$","",basename(fits))
lg("FACETS samples: %d", length(samples))
per <- list()
for (s in samples){
  af <- file.path(SPLIT, paste0(s,".tsv")); if(!file.exists(af)) next
  f <- tryCatch(readRDS(file.path(FACETS, paste0(s,"_facets_fit.rds"))), error=function(e)NULL)
  if(is.null(f) || !all(c("jointseg","out") %in% names(f))) next
  js <- as.data.table(f$jointseg); cn <- as.data.table(f$out)
  if(!all(c("tcn","lcn","seg") %in% names(cn))) next
  het <- merge(js[het==1 & keep==1, .(chrom,maploc,seg)], cn[,.(seg,tcn,lcn)], by="seg")
  het[, chrom:=as.character(chrom)]; het <- het[chrom %in% as.character(1:22)]
  a <- tryCatch(fread(af, header=FALSE, col.names=c("chrom","pos","ref","alt","refm","altm")), error=function(e)NULL)
  if(is.null(a)||nrow(a)==0) next
  a[, chrom:=as.character(chrom)]; a <- a[(refm+altm)>=10 & nchar(ref)==1 & nchar(alt)==1]
  if(nrow(a)==0) next
  a[, imb:=abs(altm/(refm+altm)-0.5)]; a[, dp:=refm+altm]
  m <- merge(het, a, by.x=c("chrom","maploc"), by.y=c("chrom","pos"))
  bal <- m[tcn==2 & lcn==1]
  loh <- m[tcn==2 & lcn==0]
  per[[s]] <- data.table(sid=sub("^T","",s),
     imb_baldiploid = if(nrow(bal)>0) weighted.mean(bal$imb,bal$dp) else NA_real_, n_bal=nrow(bal),
     imb_cnloh      = if(nrow(loh)>0) weighted.mean(loh$imb,loh$dp) else NA_real_, n_loh=nrow(loh))
}
D <- rbindlist(per); D <- D[n_bal>=20]
lg("samples with >=20 balanced-diploid het SNPs: %d", nrow(D))

cl <- fread("/path/to/data_root/CPGEA/CPGEA_web/CPGEA_Clinical_Info.tsv"); cl[, sid:=as.character(`Sequencing ID`)]
parseGS<-function(v) sapply(as.character(v),function(x){ if(grepl("\\+",x)) sum(num(strsplit(x,"\\+")[[1]])) else num(x)})
cl[, gleason:=parseGS(GS)]
pu <- fread(file.path(PROJ,"figures/supp_panels/facets_purity_ploidy_per_sample.csv")); pu[, sid:=sub("^T","",as.character(sample))]
d <- merge(D, cl[,.(sid,gleason)], by="sid"); d <- merge(d, pu[,.(sid,purity,ploidy)], by="sid")
fwrite(d, file.path(OUT,"ase_gleason_true_copyneutral.csv"))

sp<-function(a,b){dd<-d[is.finite(get(a))&is.finite(get(b))];ct<-suppressWarnings(cor.test(dd[[a]],dd[[b]],method="spearman"));sprintf("rho=%.3f P=%.2g n=%d",ct$estimate,ct$p.value,nrow(dd))}
sink(file.path(OUT,"ase_gleason_true_copyneutral_summary.txt"))
cat("TRUE copy-neutral (FACETS balanced-diploid tcn=2/lcn=1) allelic imbalance vs Gleason\n\n")
cat("  gleason ~ imb_baldiploid (regulatory):", sp("gleason","imb_baldiploid"),"\n")
cat("  gleason ~ imb_cnloh (CN-LOH)         :", sp("gleason","imb_cnloh"),"\n\n")
dd<-d[is.finite(imb_baldiploid)&is.finite(gleason)&is.finite(purity)&is.finite(ploidy)]
m<-lm(imb_baldiploid ~ gleason + purity + ploidy, dd); co<-summary(m)$coefficients
cat("=== lm: imb_baldiploid ~ gleason + purity + ploidy ===\n")
for(v in rownames(co)) cat(sprintf("  %-10s coef=%.6f P=%.3g\n", v, co[v,1], co[v,4]))
sink()
cat(readLines(file.path(OUT,"ase_gleason_true_copyneutral_summary.txt")), sep="\n")
lg("Done.")
