#!/usr/bin/env Rscript
# 58: publication version of main Fig 1 B-E (A = study design, drawn by hand; paired-comparison add-on in
# submission/fig1A_paired_schematic/). Replaces 25_pub_figs.R (B), 29_fig1d_threeset_burden.R (C), 02_Fig1C_two_genes_ASE.R (D, E).
# D is H19 (ENSG00000130600), as in v9; figure_freeze/Fig1/data/MEST_ase_data.csv is stale and not used.
# Style: fig_theme.R (178 mm, Arial 8 pt, no in-panel titles). Legend statistics -> fig1_legend_stats.txt.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
P <- "/path/to/data_root/ASE"; O <- file.path(V, "pub_figs_v9"); dir.create(O, showWarnings = FALSE)
cT <- PAL[["tumor"]]; cN <- PAL[["normal"]]
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))

pca <- fread(file.path(P, "code/PCa_ASE_genes.csv"))
note("Fig1B: %d PCa-ASE genes; tumors with significant ASE per gene: median %d, range %d-%d", nrow(pca), median(pca$n_sig_Tumor), min(pca$n_sig_Tumor), max(pca$n_sig_Tumor))
B1 <- ggplot(pca, aes(n_sig_Tumor)) + geom_histogram(binwidth = 1, fill = "#6A51A3", colour = "white", linewidth = 0.15) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) + labs(x = "Tumors with significant ASE", y = "Genes") + theme_pub()

S <- fread(file.path(P, "figure_freeze/Fig1/data/Fig1D_threeset_burden_data.csv"))
S[, set := factor(set, levels = c("All genes", "Cancer genes", "Imprinted controls"), labels = c("All", "Cancer", "Imprinted"))]
st <- S[, .(n = .N, P = wilcox.test(Tumor, Normal, paired = TRUE)$p.value, mT = median(Tumor), mN = median(Normal)), by = set]
for (i in seq_len(nrow(st))) note("Fig1C %s: n=%d genes; median |a| tumor %.3f vs normal %.3f; paired Wilcoxon P=%.2g", st$set[i], st$n[i], st$mT[i], st$mN[i], st$P[i])
L <- melt(S[, .(set, Normal, Tumor)], id.vars = "set", variable.name = "tissue", value.name = "a")
C1 <- ggplot(L, aes(set, a, fill = tissue)) + geom_boxplot(width = 0.62, linewidth = 0.3, outlier.size = 0.15, outlier.alpha = 0.3, position = position_dodge(0.7)) +
  scale_fill_manual(values = c(Normal = cN, Tumor = cT), name = NULL) + labs(x = NULL, y = "Mean ASE per gene, |a|") + theme_pub() +
  theme(legend.position = "top", legend.margin = margin(0, 0, -2, 0, "mm"))

gene <- function(gid) rbindlist(lapply(c("T", "N"), function(t) {
  fs <- list.files(P, pattern = sprintf("^%s\\d+_ase_results\\.csv$", t), full.names = TRUE)
  rbindlist(lapply(fs, function(f) fread(f)[gene_id == gid, .(a, fdr)][, sid := gsub("\\D", "", basename(f))]))[, tissue := if (t == "T") "Tumor" else "Normal"] }))[is.finite(a)]
one <- function(gid, sym) { d <- gene(gid)[, `:=`(absa = abs(a), ase = fifelse(fdr < 0.1, "ASE", "No ASE"))]
  d[, tissue := factor(tissue, levels = c("Normal", "Tumor"))]
  pr <- dcast(d, sid ~ tissue, value.var = "absa")[!is.na(Tumor) & !is.na(Normal)]
  w <- wilcox.test(pr$Tumor, pr$Normal, paired = TRUE)
  note("Fig1 %s: n tumor %d, normal %d (matched pairs %d); ASE (FDR<0.1) tumor %d, normal %d; median |a| tumor %.3f vs normal %.3f; paired Wilcoxon signed-rank P=%.2g (tumor > normal in %d of %d pairs)",
       sym, d[tissue == "Tumor", .N], d[tissue == "Normal", .N], nrow(pr), d[tissue == "Tumor", sum(ase == "ASE")], d[tissue == "Normal", sum(ase == "ASE")],
       d[tissue == "Tumor", median(absa)], d[tissue == "Normal", median(absa)], w$p.value, sum(pr$Tumor > pr$Normal), nrow(pr))
  set.seed(58)
  ggplot(d, aes(tissue, absa)) + geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(aes(colour = ase), width = 0.15, height = 0, size = 0.6) +
    scale_colour_manual(values = c(ASE = "#D6604D", `No ASE` = "#92C5DE"), name = NULL) +
    scale_y_continuous(limits = c(0, 0.5)) + labs(x = NULL, y = bquote(italic(.(sym))~"|a|")) + theme_pub() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -2, 0, "mm")) }
D1 <- one("ENSG00000130600", "H19"); E1 <- one("ENSG00000073282", "TP63")

g <- (B1 | C1 | D1 | E1) + plot_layout(widths = c(1.1, 1.5, 0.8, 0.8)) + plot_annotation(tag_levels = list(c("B", "C", "D", "E"))) & tags
save_panel(g, file.path(O, "Fig1_BE.pdf"), 178, 58)
writeLines(LS, file.path(O, "fig1_legend_stats.txt")); cat("DONE\n")
