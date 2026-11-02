#!/usr/bin/env Rscript
# 58b (2026-10-13): v10 Fig 1 C-F in one 178 x 62 mm row (A, B = hand-drawn schematics, assembled separately).
# Changes from 58: P values written as m x 10^e on the panels (no R-style 2.2e-16), ASE status in greys so red/blue mean tumor/normal only,
# panel letters C-F as in the v10 legend.
# 58: publication version of main Fig 1 B-E (A = study design, drawn by hand; paired-comparison add-on in
# submission/fig1A_paired_schematic/). Replaces 25_pub_figs.R (B), 29_fig1d_threeset_burden.R (C), 02_Fig1C_two_genes_ASE.R (D, E).
# D is H19 (ENSG00000130600), as in v9; figure_freeze/Fig1/data/MEST_ase_data.csv is stale and not used.
# Style: fig_theme.R (178 mm, Arial 8 pt, no in-panel titles). Legend statistics -> fig1_legend_stats.txt.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
P <- "/path/to/data_root/ASE"; O <- file.path(V, "pub_figs_v10"); dir.create(O, showWarnings = FALSE)
cT <- PAL[["tumor"]]; cN <- PAL[["normal"]]
pl <- function(p) if (p == 0) "italic(P) < 2.2 %*% 10^-16" else if (p < 1e-3) { e <- floor(log10(p)); sprintf("italic(P) == %.1f %%*%% 10^%d", p / 10^e, e) } else sprintf("italic(P) == %s", formatC(p, digits = 2, format = "fg"))
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))

pca <- fread(file.path(P, "code/PCa_ASE_genes.csv"))
note("Fig1B: %d PCa-ASE genes; tumors with significant ASE per gene: median %d, range %d-%d", nrow(pca), median(pca$n_sig_Tumor), min(pca$n_sig_Tumor), max(pca$n_sig_Tumor))
B1 <- ggplot(pca, aes(n_sig_Tumor)) + geom_histogram(binwidth = 1, fill = "#6A51A3", colour = "white", linewidth = 0.15) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) + labs(x = "Tumors with significant ASE", y = "PCa-ASE genes") + theme_pub()

S <- fread(file.path(P, "figure_freeze/Fig1/data/Fig1D_threeset_burden_data.csv"))
S[, set := factor(set, levels = c("All genes", "Cancer genes", "Imprinted controls"), labels = c("All", "Cancer", "Imprinted"))]
st <- S[, .(n = .N, P = wilcox.test(Tumor, Normal, paired = TRUE)$p.value, mT = median(Tumor), mN = median(Normal)), by = set]
for (i in seq_len(nrow(st))) note("Fig1C %s: n=%d genes; median |a| tumor %.3f vs normal %.3f; paired Wilcoxon P=%.2g", st$set[i], st$n[i], st$mT[i], st$mN[i], st$P[i])
L <- melt(S[, .(set, Normal, Tumor)], id.vars = "set", variable.name = "tissue", value.name = "a")
C1 <- ggplot(L, aes(set, a, fill = tissue)) + geom_boxplot(width = 0.62, linewidth = 0.3, outlier.size = 0.15, outlier.alpha = 0.3, position = position_dodge(0.7)) +
  geom_text(data = st[, lab := fifelse(P < 0.001, "***", fifelse(P < 0.01, "**", fifelse(P < 0.05, "*", "ns")))], aes(set, fifelse(lab == "ns", 0.55, 0.53), label = lab),
            inherit.aes = FALSE, size = BASE / .pt, family = FONT) +
  scale_fill_manual(values = c(Normal = cN, Tumor = cT), name = NULL) + scale_y_continuous(limits = c(0, 0.58), breaks = seq(0, 0.5, 0.1)) +
  labs(x = NULL, y = expression("Mean ASE per gene, |"*italic(a)*"|")) + theme_pub() +
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
    annotate("text", x = 1.5, y = 0.56, label = pl(w$p.value), parse = TRUE, size = (BASE - 1) / .pt, family = FONT) +
    scale_colour_manual(values = c(ASE = "black", `No ASE` = "grey70"), labels = c(ASE = "ASE (FDR < 0.1)", `No ASE` = "No ASE"), name = NULL) +
    scale_y_continuous(limits = c(0, 0.58), breaks = seq(0, 0.5, 0.1)) + labs(x = NULL, y = bquote(italic(.(sym))*", |"*italic(a)*"|")) + theme_pub() +
    theme(legend.position = "top", legend.margin = margin(0, 0, -2, 0, "mm")) }
D1 <- one("ENSG00000130600", "H19"); E1 <- one("ENSG00000073282", "TP63")

g <- ((B1 | C1 | D1 | E1) + plot_layout(widths = c(1.1, 1.5, 0.8, 0.8), guides = "keep") + plot_annotation(tag_levels = list(c("C", "D", "E", "F"))) & tags)
g[[4]] <- g[[4]] + theme(legend.position = "none")
save_panel(g, file.path(O, "Fig1_CF.pdf"), 178, 62)
writeLines(LS, file.path(O, "fig1_legend_stats.txt")); cat("DONE\n")
