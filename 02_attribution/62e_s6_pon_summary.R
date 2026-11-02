#!/usr/bin/env Rscript
# 62e: S6 data on option A (62d: PoN-filtered, ANNOVAR, ASE cohort, no gnomAD filter). Consequence counts once per sample-gene-consequence;
# Venn of missense-carrying (nonsynonymous SNV) vs SCNA-associated PCa-ASE genes (gene symbols). Reported with and without T438/T502
# (279,045 and 195,210 PoN-filtered calls vs ~5-30k typical: likely tumor-normal or normal-sample problems).
suppressPackageStartupMessages(library(data.table))
V <- "/path/to/data_root/ASE/code/v8_round2"
X <- fread(file.path(V, "supp_data/S6_pon_exonic_variants.tsv"))[, gene := sub("[,;].*", "", gene)]
pca <- fread("/path/to/data_root/ASE/code/PCa_ASE_genes.csv")
scna <- fread("/path/to/data_root/ASE/code/fig2_results/pca_ase_cnv_correlations.csv")[is.finite(fdr) & fdr < 0.1]
lab <- c(nonsynonymous_SNV = "Missense", stopgain = "Stop gained", frameshift_deletion = "Frameshift", frameshift_insertion = "Frameshift",
         nonframeshift_deletion = "Inframe deletion", nonframeshift_insertion = "Inframe insertion", nonframeshift_substitution = "Inframe substitution",
         startloss = "Start lost", stoploss = "Stop lost")
X[, cl := fifelse(func %like% "splicing" & !(func %like% "exonic"), "Splice site", lab[exonic_func])]
pc <- X[, .N, by = sample][order(-N)]; cat("exonic/splicing PoN-filtered calls per tumor: median", median(pc$N), "; top:", paste(head(pc$sample, 4), head(pc$N, 4), collapse = ", "), "\n")
for (excl in list(character(0), c("T438", "T502"))) {
  Y <- X[!sample %in% excl & gene %in% pca$gene_symbol & !is.na(cl)]
  cnt <- unique(Y[, .(sample, gene, cl)])[, .N, by = cl][order(-N)]
  mg <- unique(Y[cl == "Missense", gene]); ov <- sum(scna$gene_symbol %in% mg)
  tag <- if (length(excl)) "excl_T438_T502" else "all"
  cat(sprintf("[%s] tumors %d ; consequence counts: %s\n  missense-carrying PCa-ASE genes %d ; SCNA-associated %d ; overlap %d\n", tag, uniqueN(X[!sample %in% excl, sample]),
      cnt[, paste(cl, N, collapse = ", ")], length(mg), nrow(scna), ov))
  fwrite(cnt, file.path(V, sprintf("supp_data/S6A_consequence_counts_PoN_%s.csv", tag)))
  fwrite(data.table(set = c("SCNA", "missense", "overlap"), n = c(nrow(scna), length(mg), ov)), file.path(V, sprintf("supp_data/S6C_venn_PoN_%s.csv", tag)))
}
