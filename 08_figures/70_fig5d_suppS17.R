#!/usr/bin/env Rscript
# 70: publication versions of (i) main Fig 5D, FOXM1 -> KIF18B expression in tumors (model and inputs as code/fig4b_scatter.R:
#     lm(log2 KIF18B ~ FOXM1 z + FACETS purity + 6-gene proliferation score)), and (ii) Supplemental Fig S17 (v10), event-level source
#     attribution and confidence tier (saved counts from 02_event_attribution_domseg.R, as plotted by 03_S10_replot_domseg.R).
# Style: fig_theme.R. Legend numbers -> fig5d_suppS17_legend_stats.txt.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
O <- file.path(V, "pub_figs_v9")
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
TS <- BASE / .pt; tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))

fpkm <- fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls")
cols <- grep("^T\\d+_WTS$", names(fpkm), value = TRUE)
ex <- function(sym) data.table(sample_id = sub("_WTS$", "", sub("^T", "", cols)), v = as.numeric(fpkm[gene_name == sym][1, ..cols]))
PROLIF <- c("MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "BUB1")
e <- Reduce(function(a, b) merge(a, b, by = "sample_id"), lapply(c("KIF18B", "FOXM1", PROLIF), function(s) setnames(ex(s), "v", s)))
e[, prolif := rowMeans(scale(log2(as.matrix(e[, ..PROLIF]) + 0.01)))]
pur <- fread("/path/to/data_root/variants/batch_43/facets/Table_S1_PRAD.txt")[, .(sample_id = as.character(SAMPLE_ID), purity = as.numeric(purity))]
coh <- gsub("\\D", "", list.files("/path/to/data_root/ASE", pattern = "^T\\d+_ase_results\\.csv$")); stopifnot(length(coh) == 118L)
pur <- pur[sample_id %in% coh]; note("Fig5D: purity imputed (cohort median %.3f) for %s", median(pur$purity, na.rm = TRUE), paste(pur[is.na(purity), sample_id], collapse = ","))
pur[is.na(purity), purity := median(pur$purity, na.rm = TRUE)]
d <- merge(e, pur, by = "sample_id")[, `:=`(y = log2(KIF18B + 0.01), FOXM1z = as.numeric(scale(log2(FOXM1 + 0.01))))]
co <- summary(lm(y ~ FOXM1z + purity + prolif, d))$coefficients
ct <- cor.test(d$FOXM1z, d$y, method = "s", exact = FALSE)
note("Fig5D: n=%d tumors; adjusted beta (per SD FOXM1) %.2f (SE %.3f), P=%.2g; unadjusted Spearman rho %.2f, P=%.2g",
     nrow(d), co["FOXM1z", 1], co["FOXM1z", 2], co["FOXM1z", 4], ct$estimate, ct$p.value)
D5 <- ggplot(d, aes(FOXM1z, y)) + geom_point(size = 0.9, colour = PAL[["dark"]], alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, colour = PAL[["tumor"]], fill = PAL[["tumor"]], alpha = 0.15, linewidth = 0.5) +
  labs(x = expression(italic("FOXM1")*", log"[2]*" FPKM, z"), y = expression(italic("KIF18B")*", log"[2]*" FPKM")) + theme_pub()
save_panel(D5, file.path(O, "Fig5D_FOXM1_KIF18B.pdf"), 70, 60)

R <- file.path(V, "event_attribution_domseg")
ep <- fread(file.path(R, "event_level_source_proportions.csv"))
det <- fread(file.path(R, "event_level_attribution.csv")); cp <- det[, .(n_events = .N), by = source_confidence][, pct := 100 * n_events / sum(n_events)]
lab <- c(SCNA_deletion = "Deletion or LOH", copy_neutral_regulatory = "Balanced copy-neutral", unresolved = "Unresolved",
         SCNA_amplification = "Imbalanced amplification", CN_LOH = "Copy-neutral LOH", balanced_gain = "Balanced gain")
stopifnot(setequal(ep$putative_source, names(lab)))
ep[, src := factor(lab[putative_source], levels = lab[ep[order(-n_events), putative_source]])]
note("S17A: %d events; %s", sum(ep$n_events), ep[order(-n_events), paste(sprintf("%s %d (%.1f%%)", lab[putative_source], n_events, pct_events), collapse = "; ")])
note("S17B: %s", cp[order(match(source_confidence, c("HIGH", "MODERATE", "LOW"))), paste(sprintf("%s %d (%.1f%%)", source_confidence, n_events, pct), collapse = "; ")])
fillA <- c(`Deletion or LOH` = "#2166AC", `Imbalanced amplification` = PAL[["tumor"]], `Copy-neutral LOH` = "#4393C3",
           `Balanced copy-neutral` = "#F4A582", `Balanced gain` = "#FDDBC7", Unresolved = "grey70")
A17 <- ggplot(ep, aes(src, n_events, fill = src)) + geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_events)), vjust = -0.3, size = TS, family = FONT) +
  scale_fill_manual(values = fillA, guide = "none") + scale_y_continuous(expand = expansion(c(0, 0.1))) +
  labs(x = NULL, y = "ASE events") + theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
cp[, tier := factor(tools::toTitleCase(tolower(source_confidence)), levels = c("High", "Moderate", "Low"))]
B17 <- ggplot(cp, aes(tier, n_events, fill = tier)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.3, size = TS, family = FONT) +
  scale_fill_manual(values = c(High = "#1A9850", Moderate = "#FEE08B", Low = "grey75"), guide = "none") +
  scale_y_continuous(expand = expansion(c(0, 0.1))) + labs(x = "Source confidence", y = "ASE events") + theme_pub()
save_panel((A17 | B17) + plot_layout(widths = c(2, 1)) + plot_annotation(tag_levels = "A") & tags, file.path(O, "SuppFig_S17_event.pdf"), 150, 70)
writeLines(LS, file.path(O, "fig5d_suppS17_legend_stats.txt")); cat("DONE\n")
