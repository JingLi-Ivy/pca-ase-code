#!/usr/bin/env Rscript
suppressMessages(library(data.table)); setDTthreads(8)
BASE <- "/path/to/data_root/ASE"; OUT <- file.path(BASE,"code","paired_diff_ase")
allres <- fread(file.path(OUT,"gw_depthmatched_patient_gene.csv.gz"))
cat("rows:", nrow(allres), " depths:", paste(sort(unique(allres$D)),collapse=","), "\n")
pca <- fread(file.path(BASE,"code","PCa_ASE_genes.csv")); recur <- unique(sub("\\..*","",pca$gene_id))
map <- fread(file.path(BASE,"code","ensembl_to_symbol.csv")); setnames(map,1:2,c("gid","symbol"))
cand <- c("KIF18B","NEK2","TOP2A","CELF2","EVI5","ITPR1","DRD1","CNN1","GLB1L3")
summ <- list()
for (dd in sort(unique(allres$D))){
  sub <- allres[D==dd]
  r <- sub[, { n<-.N; if(n<10) NULL else { wt<-suppressWarnings(wilcox.test(aT,aN,paired=TRUE))
      .(n_pair=n, med_dT=median(aT), med_dN=median(aN), delta=median(aT-aN), p=wt$p.value) } }, by=gid]
  r <- r[!is.na(p)]
  r[, fdr := p.adjust(p,"BH")]
  r[, dir := fifelse(delta>0,"tumor_up", fifelse(delta<0,"normal_up","no_change"))]
  r[, in_recurrence := gid %in% recur]
  r <- merge(r, map, by="gid", all.x=TRUE); setorder(r, p); r[, rank_by_p := .I]
  fwrite(r, file.path(OUT, sprintf("gw_depthmatched_pergene_D%d.csv", dd)))
  cat(sprintf("\n===== D=%d : tested %d | FDR<0.1 tumor_up %d | normal_up %d | recur testable %d confirmed %d discordant %d\n",
    dd, nrow(r), r[fdr<0.1&delta>0,.N], r[fdr<0.1&delta<0,.N],
    r[in_recurrence==TRUE,.N], r[in_recurrence&fdr<0.1&delta>0,.N], r[in_recurrence&fdr<0.1&delta<0,.N]))
  print(r[symbol %in% cand, .(symbol,n_pair,medT=round(med_dT,3),medN=round(med_dN,3),
        delta=round(delta,4),p=signif(p,3),fdr=signif(fdr,3),dir)][order(-delta)])
  summ[[as.character(dd)]] <- data.table(D=dd, genes_tested=nrow(r),
    cat_tumor_up=r[fdr<0.1&delta>0,.N], cat_normal_up=r[fdr<0.1&delta<0,.N],
    recur_testable=r[in_recurrence==TRUE,.N], recur_confirmed=r[in_recurrence&fdr<0.1&delta>0,.N],
    recur_discordant=r[in_recurrence&fdr<0.1&delta<0,.N],
    pct_confirmed=round(100*r[in_recurrence&fdr<0.1&delta>0,.N]/max(r[in_recurrence==TRUE,.N],1),1))
}
S <- rbindlist(summ); fwrite(S, file.path(OUT,"gw_depthmatched_summary.csv"))
cat("\n===== SUMMARY =====\n"); print(S); cat("DONE\n")
