#!/usr/bin/env Rscript
# 69: Supplemental Fig S18 (v10) - GO BP enrichment of PCa-ASE genes, redrawn from 68 (topGO classic Fisher, nodeSize 5, background = genes with
# heterozygous-SNP coverage; BH over ALL tested terms). Shows every term with FDR < 0.10, ordered by P; colour = -log10 FDR, size = PCa-ASE genes.
# Replaces Fig1F_GO_enrichment (BH over top-100 terms only).
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R")); O <- file.path(V, "pub_figs_v9")
G <- fread(file.path(V, "supp_data/GO_BP_topGO_classic_allterms_FDR.csv"))[FDR_all < 0.10][order(p)]
G[, Term := factor(Term, levels = rev(Term))]
hl <- c("GO:0016192", "GO:0070972", "GO:0048193")
p <- ggplot(G, aes(ratio, Term)) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_point(aes(size = Significant, colour = -log10(FDR_all))) +
  scale_colour_gradient(low = "#9ECAE1", high = PAL[["tumor"]], name = expression(-log[10]~FDR)) +
  scale_size_area(max_size = 3.5, breaks = c(20, 100, 500), name = "PCa-ASE genes") +
  labs(x = "Enrichment ratio (observed / expected)", y = NULL) + theme_pub() +
  theme(axis.text.y = element_text(face = ifelse(rev(G$GO.ID) %in% hl, "bold", "plain")))
save_panel(p, file.path(O, "SuppFig_S18_GO.pdf"), 150, 150)
L <- c(sprintf("terms tested %d; FDR<0.10 %d; FDR<0.05 %d", nrow(fread(file.path(V, "supp_data/GO_BP_topGO_classic_allterms_FDR.csv"))), nrow(G), G[FDR_all < 0.05, .N]),
       G[GO.ID %in% hl, sprintf("%s %s: %d/%.1f genes, ratio %.2f, P %.2g, FDR %.2g, rank %d", GO.ID, Term, Significant, Expected, ratio, p, FDR_all, rank_p)])
writeLines(L, file.path(O, "SuppFig_S18_legend_stats.txt")); cat(L, sep = "\n"); cat("DONE\n")
