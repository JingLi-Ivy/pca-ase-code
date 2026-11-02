#!/usr/bin/env Rscript
# Robustness for ASE(mean_imb)-clinical: is Gleason/BCR link real vs purity/aneuploidy/grade?
suppressPackageStartupMessages({ library(data.table); library(survival) })
PROJ <- "/path/to/data_root/ASE"; OUT <- file.path(PROJ,"figures/ase_clinical")
lg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), sprintf(...)))
num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]","",as.character(x))))

ev <- fread(file.path(PROJ,"code/event_attribution_results/event_level_attribution.csv"))
ev[, sid:=sub("^T","",as.character(Sample))]; ev[, gid:=sub("\\..*","",gene_id_clean)]
ev[, cn:=ifelse(cnv_type=="NEUTRAL","neutral","scna")]
cnmap <- unique(ev[,.(sid,gid,cn)])
ase <- fread(file.path(PROJ,"all_T_samples_SNPs_with_ASE_annotation.csv"))
ase[, sid:=sub("^T","",as.character(sample_id))]; ase[, gid:=sub("\\..*","",gene_id)]
ase <- ase[(ref.matches+alt.matches)>=10 & ref_ratio>0 & ref_ratio<1]
ase[, imb:=abs(ref_ratio-0.5)]; ase[, dp:=ref.matches+alt.matches]
ase <- merge(ase, cnmap, by=c("sid","gid"), all.x=TRUE); ase[is.na(cn), cn:="neutral_or_untested"]
b <- ase[, .(imb_all=weighted.mean(imb,dp)), by=sid]
bn <- ase[cn!="scna", .(imb_neutral=weighted.mean(imb,dp)), by=sid]
bs <- ase[cn=="scna", .(imb_scna=weighted.mean(imb,dp)), by=sid]

cl <- fread("/path/to/data_root/CPGEA/CPGEA_web/CPGEA_Clinical_Info.tsv"); cl[, sid:=as.character(`Sequencing ID`)]
parseGS <- function(v) sapply(as.character(v), function(s){ if(grepl("\\+",s)) sum(num(strsplit(s,"\\+")[[1]])) else num(s)})
cl[, `:=`(gleason=parseGS(GS), psa=num(PSA), logpsa=log2(num(PSA)+1),
          bcr=fifelse(toupper(as.character(BCR))=="YES",1L,fifelse(toupper(as.character(BCR))=="NO",0L,NA_integer_)),
          months=num(`Month to surgery/BCR`))]
pu <- fread(file.path(PROJ,"figures/supp_panels/facets_purity_ploidy_per_sample.csv")); pu[, sid:=sub("^T","",as.character(sample))]
d <- Reduce(function(x,y) merge(x,y,by="sid",all.x=TRUE), list(b,bn,bs,
      cl[,.(sid,gleason,psa,logpsa,bcr,months)], pu[,.(sid,purity,ploidy)]))
lg("samples: %d", nrow(d))
fwrite(d, file.path(OUT,"ase_clinical_robust_data.csv"))

sink(file.path(OUT,"ase_clinical_robust_summary.txt"))
sp <- function(a,b){ dd<-d[is.finite(get(a))&is.finite(get(b))]; ct<-suppressWarnings(cor.test(dd[[a]],dd[[b]],method="spearman")); sprintf("rho=%.3f P=%.2g n=%d",ct$estimate,ct$p.value,nrow(dd)) }
cat("=== Gleason vs ASE imbalance (decomposed) ===\n")
cat("  gleason ~ imb_all      :", sp("gleason","imb_all"),"\n")
cat("  gleason ~ imb_neutral  :", sp("gleason","imb_neutral"),"  <- copy-neutral(regulatory) component\n")
cat("  gleason ~ imb_scna     :", sp("gleason","imb_scna"),"  <- SCNA-region component\n")
cat("  gleason ~ purity       :", sp("gleason","purity"),"\n")
cat("  gleason ~ ploidy       :", sp("gleason","ploidy"),"\n")
cat("  imb_all ~ purity       :", sp("imb_all","purity"),"\n")

cat("\n=== lm: imb_all ~ gleason + purity + ploidy (is grade effect robust?) ===\n")
m <- lm(imb_all ~ gleason + purity + ploidy, data=d); co<-summary(m)$coefficients
for(v in rownames(co)) cat(sprintf("  %-12s coef=%.5f P=%.3g\n", v, co[v,1], co[v,4]))
cat("\n=== lm: imb_neutral ~ gleason + purity + ploidy (regulatory component) ===\n")
m2 <- lm(imb_neutral ~ gleason + purity + ploidy, data=d); co2<-summary(m2)$coefficients
for(v in rownames(co2)) cat(sprintf("  %-12s coef=%.5f P=%.3g\n", v, co2[v,1], co2[v,4]))

cat("\n=== BCR Cox: univariate vs multivariable (mean_imb per-SD) ===\n")
ds <- d[is.finite(months)&months>0&!is.na(bcr)&is.finite(imb_all)&is.finite(gleason)&is.finite(logpsa)]
ds[, z:=as.numeric(scale(imb_all))]
u <- summary(coxph(Surv(months,bcr)~z, ds))
mv<- summary(coxph(Surv(months,bcr)~z+gleason+logpsa, ds))
cat(sprintf("  univariate       : HR=%.3f (%.2f-%.2f) P=%.3g\n", u$conf.int[1,1],u$conf.int[1,3],u$conf.int[1,4],u$coefficients[1,5]))
cat(sprintf("  +Gleason+PSA(adj) : imb HR=%.3f (%.2f-%.2f) P=%.3g | gleason P=%.3g | logpsa P=%.3g (n=%d,ev=%d)\n",
    mv$conf.int["z",1],mv$conf.int["z",3],mv$conf.int["z",4],mv$coefficients["z",5],
    mv$coefficients["gleason",5],mv$coefficients["logpsa",5],nrow(ds),sum(ds$bcr)))
ds2 <- ds[is.finite(imb_neutral)]; ds2[, zn:=as.numeric(scale(imb_neutral))]
mv2<- summary(coxph(Surv(months,bcr)~zn+gleason+logpsa, ds2))
cat(sprintf("  copy-neutral imb +Gleason+PSA: HR=%.3f (%.2f-%.2f) P=%.3g\n",
    mv2$conf.int["zn",1],mv2$conf.int["zn",3],mv2$conf.int["zn",4],mv2$coefficients["zn",5]))
sink()
cat(readLines(file.path(OUT,"ase_clinical_robust_summary.txt")), sep="\n")
lg("Done.")
