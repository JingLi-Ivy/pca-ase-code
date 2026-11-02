#!/usr/bin/env Rscript
## Genome-wide paired differential ASE with depth control (depth-matched downsampling,
## equal-coverage restriction, threshold sensitivity)
##
## Design:
##  Each patient's T/N pair is processed separately; only heterozygous sites shared by both tissues (same chr:pos) are kept,
##  then limited by a per-site depth threshold D, and every site in both tissues is downsampled without replacement to exactly D reads
##  (rhyper), removing both depth and site-composition confounding. Gene-level statistic = mean site |ratio-0.5|.
##  Results are averaged over R repeats to reduce sampling noise. Finally a paired Wilcoxon test per gene, then BH.
suppressMessages({library(data.table); library(parallel)})
setDTthreads(2)
BASE <- "/path/to/data_root/ASE"
OUT  <- file.path(BASE,"code","paired_diff_ase"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
NREP <- 5L
DEPTHS <- c(8L, 12L, 20L)
set.seed(42)

gb <- fread(file.path(BASE,"gtf_gene_body_geneid.bed"), header=FALSE,
            col.names=c("chr","start","end","gid"))
gb[, chr := as.character(chr)]
gb <- gb[!is.na(start) & !is.na(end) & end>=start]
setkey(gb, chr, start, end)
cat(sprintf("gene intervals: %d\n", nrow(gb)))

pids <- sub("^T","", sub("_allele.*","", basename(
  Sys.glob(file.path(BASE,"T*_allele_specific_counts.pipeline.nonzero.txt")))))
pids <- pids[file.exists(file.path(BASE, sprintf("N%s_allele_specific_counts.pipeline.nonzero.txt", pids)))]
cat(sprintf("paired patients: %d\n", length(pids)))

rd <- function(f){
  x <- fread(f, header=FALSE, showProgress=FALSE, select=1:7,
             col.names=c("chr","pos","ref","alt","rs","nref","nalt"))
  x[, chr := as.character(chr)]
  x[, dp := nref + nalt]
  x[dp > 0, .(chr,pos,nref,nalt,dp)]
}

one <- function(pid){
  fT <- file.path(BASE, sprintf("T%s_allele_specific_counts.pipeline.nonzero.txt", pid))
  fN <- file.path(BASE, sprintf("N%s_allele_specific_counts.pipeline.nonzero.txt", pid))
  tt <- tryCatch(rd(fT), error=function(e) NULL); nn <- tryCatch(rd(fN), error=function(e) NULL)
  if (is.null(tt) || is.null(nn)) return(NULL)
  m <- merge(tt, nn, by=c("chr","pos"), suffixes=c("_T","_N"))
  if (!nrow(m)) return(NULL)
  m[, `:=`(start=pos, end=pos)]
  setkey(m, chr, start, end)
  ov <- foverlaps(m, gb, type="within", nomatch=0L)
  if (!nrow(ov)) return(NULL)
  res <- rbindlist(lapply(DEPTHS, function(D){
    s <- ov[dp_T >= D & dp_N >= D]
    if (!nrow(s)) return(NULL)
    accT <- numeric(nrow(s)); accN <- numeric(nrow(s))
    for (r in seq_len(NREP)){
      kT <- rhyper(nrow(s), s$nref_T, s$nalt_T, D)
      kN <- rhyper(nrow(s), s$nref_N, s$nalt_N, D)
      accT <- accT + abs(kT/D - 0.5); accN <- accN + abs(kN/D - 0.5)
    }
    s[, `:=`(aT = accT/NREP, aN = accN/NREP)]
    s[, .(n_site=.N, aT=mean(aT), aN=mean(aN)), by=.(gid)][, `:=`(pid=pid, D=D)]
  }))
  res
}

t0 <- Sys.time()
allres <- rbindlist(mclapply(pids, function(p) tryCatch(one(p), error=function(e){
  cat("FAIL", p, conditionMessage(e), "\n"); NULL}), mc.cores=16), fill=TRUE)
cat(sprintf("per-patient stage done in %.1f min; rows=%d\n",
            as.numeric(difftime(Sys.time(), t0, units="mins")), nrow(allres)))
fwrite(allres, file.path(OUT,"gw_depthmatched_patient_gene.csv.gz"))

gtest <- function(d){
  d[, {
    n <- .N
    if (n < 10) NULL else {
      wt <- suppressWarnings(wilcox.test(aT, aN, paired=TRUE))
      .(n_pair=n, med_dT=median(aT), med_dN=median(aN),
        delta=median(aT-aN), p=wt$p.value)
    }
  }, by=gid]
}
pca <- fread(file.path(BASE,"code","PCa_ASE_genes.csv"))
recur <- unique(sub("\\..*","", pca$gene_id))
map <- fread(file.path(BASE,"code","ensembl_to_symbol.csv")); setnames(map,1:2,c("gid","symbol"))

summ <- list()
for (D in DEPTHS){
  r <- gtest(allres[D==..D])
  if (is.null(r) || !nrow(r)) next
  r[, fdr := p.adjust(p,"BH")]
  r[, dir := fifelse(delta>0,"tumor_up", fifelse(delta<0,"normal_up","no_change"))]
  r[, in_recurrence := gid %in% recur]
  r <- merge(r, map, by="gid", all.x=TRUE)
  setorder(r, p); r[, rank_by_p := .I]
  fwrite(r, file.path(OUT, sprintf("gw_depthmatched_pergene_D%d.csv", D)))
  summ[[as.character(D)]] <- data.table(
    D=D, genes_tested=nrow(r),
    catalogue_tumor_up=r[fdr<0.1 & delta>0, .N],
    catalogue_normal_up=r[fdr<0.1 & delta<0, .N],
    recur_testable=r[in_recurrence==TRUE, .N],
    recur_confirmed=r[in_recurrence & fdr<0.1 & delta>0, .N],
    recur_discordant=r[in_recurrence & fdr<0.1 & delta<0, .N])
  cat(sprintf("D=%d: tested %d | tumor_up %d | normal_up %d | recur testable %d confirmed %d discordant %d\n",
      D, nrow(r), r[fdr<0.1 & delta>0,.N], r[fdr<0.1 & delta<0,.N],
      r[in_recurrence==TRUE,.N], r[in_recurrence & fdr<0.1 & delta>0,.N],
      r[in_recurrence & fdr<0.1 & delta<0,.N]))
  cand <- c("KIF18B","NEK2","TOP2A","CELF2","EVI5","ITPR1","DRD1","CNN1","GLB1L3")
  cat("  9 candidates:\n")
  print(r[symbol %in% cand, .(symbol,n_pair,medT=round(med_dT,3),medN=round(med_dN,3),
                              delta=round(delta,4),p=signif(p,3),fdr=signif(fdr,3),dir)][order(-delta)])
}
S <- rbindlist(summ); fwrite(S, file.path(OUT,"gw_depthmatched_summary.csv")); print(S)
cat("DONE\n")
