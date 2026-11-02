#!/usr/bin/env Rscript
# Paired differential-ASE test.
# Replaces recurrence-only definition with a FORMAL paired tumor-vs-normal test of
# allelic-imbalance magnitude. Unit = gene; paired across the 118 matched patients.
#  - per-gene paired Wilcoxon signed-rank on |a| (tumor vs normal), FDR (BH)
#  - global LMM |a| ~ tissue + (1|patient) + (1|gene) for overall shift
#  - threshold sensitivity (min pairs) ; overlap/discordance vs recurrence 2378
suppressMessages({library(data.table)})
OUT <- "/path/to/data_root/ASE/code/paired_diff_ase"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
log <- function(...) cat(sprintf(...),"\n")

load("/path/to/data_root/ASE/ase_data.RData"); setDT(ase_data)
ase_data[, gid := sub("\\..*","",gene_id)]
ase_data[, pid := sub("^[TN]","", sub("_.*","", sample_id))]
ase_data[, absA := abs(a)]
d <- ase_data[, .(absA=max(absA,na.rm=TRUE), sig=any(fdr<0.1,na.rm=TRUE)),
              by=.(pid,gid,sample_type)]
d <- d[is.finite(absA)]
w <- dcast(d, pid+gid ~ sample_type, value.var="absA")
setnames(w, c("Tumor","Normal"), c("aT","aN"), skip_absent=TRUE)
w <- w[is.finite(aT) & is.finite(aN)]
log("paired patient-gene observations (evaluable both tissues): %d over %d genes",
    nrow(w), length(unique(w$gid)))

gtest <- function(sub){
  n <- nrow(sub); if (n < 10) return(NULL)
  wt <- suppressWarnings(wilcox.test(sub$aT, sub$aN, paired=TRUE))
  data.table(n_pair=n, med_dT=median(sub$aT), med_dN=median(sub$aN),
             delta=median(sub$aT - sub$aN), p=wt$p.value)
}
res <- w[, gtest(.SD), by=gid]
res <- res[!is.na(p)]
res[, fdr := p.adjust(p, method="BH")]
res[, dir := fifelse(delta>0,"tumor_up","normal_up")]
log("genes tested (>=10 pairs): %d", nrow(res))

pca <- fread("/path/to/data_root/ASE/code/PCa_ASE_genes.csv")
recur <- unique(sub("\\..*","",pca$gene_id))
res[, in_recurrence := gid %in% recur]

sig <- res[fdr<0.1 & delta>0]
log("=== PRIMARY differential-ASE catalogue (FDR<0.1 & tumor>normal): %d genes ===", nrow(sig))
log("  of which in recurrence-2378: %d ; novel (not in 2378): %d", sum(sig$in_recurrence), sum(!sig$in_recurrence))
log("recurrence-2378 genes testable here: %d", res[in_recurrence==TRUE, .N])
log("  recurrence genes confirmed by paired test (FDR<0.1 tumor_up): %d (%.1f%% of testable)",
    res[in_recurrence & fdr<0.1 & delta>0, .N],
    100*res[in_recurrence & fdr<0.1 & delta>0, .N]/max(res[in_recurrence==TRUE,.N],1))
log("  recurrence genes DISCORDANT (paired normal_up @FDR<0.1): %d",
    res[in_recurrence & fdr<0.1 & delta<0, .N])

for (mn in c(10,15,20,30)){
  rr <- res[n_pair>=mn]; rr[, fdr2:=p.adjust(p,"BH")]
  log("  sensitivity min_pairs>=%d: tested %d, FDR<0.1 tumor_up %d", mn, nrow(rr), rr[fdr2<0.1 & delta>0,.N])
}

if (requireNamespace("lme4", quietly=TRUE)){
  dd <- copy(d); dd[, tissue := factor(sample_type, levels=c("Normal","Tumor"))]
  dd <- dd[gid %in% recur]
  fit <- tryCatch(lme4::lmer(absA ~ tissue + (1|pid) + (1|gid), data=dd),
                  error=function(e){log("lmer fail: %s",conditionMessage(e));NULL})
  if(!is.null(fit)){
    ce <- summary(fit)$coefficients
    log("global LMM |a|~tissue+(1|patient)+(1|gene): tumor beta=%.4f (SE %.4f, t=%.2f)",
        ce["tissueTumor",1], ce["tissueTumor",2], ce["tissueTumor",3])
    fwrite(data.table(term=rownames(ce), est=ce[,1], se=ce[,2], t=ce[,3]),
           file.path(OUT,"global_lmm_coef.csv"))
  }
}

fwrite(res[order(fdr)], file.path(OUT,"paired_diff_ase_pergene.csv"))
fwrite(sig[order(fdr)], file.path(OUT,"differential_ase_catalogue_FDR10.csv"))
summ <- data.table(
  metric=c("paired_obs","genes_tested","primary_catalogue_FDR10_tumorUp",
           "catalogue_in_recurrence","catalogue_novel","recurrence_testable",
           "recurrence_confirmed_FDR10","recurrence_discordant_FDR10"),
  value=c(nrow(w), nrow(res), nrow(sig), sum(sig$in_recurrence), sum(!sig$in_recurrence),
          res[in_recurrence==TRUE,.N], res[in_recurrence & fdr<0.1 & delta>0,.N],
          res[in_recurrence & fdr<0.1 & delta<0,.N]))
fwrite(summ, file.path(OUT,"paired_diff_ase_summary.csv")); print(summ)
log("DONE -> %s", OUT)
