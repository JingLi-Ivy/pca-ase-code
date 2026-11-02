#!/usr/bin/env Rscript
# 64: Supplemental Fig S6 (somatic mutations in PCa-ASE genes; former v9 Fig 3), on one consistent definition:
# option A (chosen 2026-10-08): panel-of-normals-filtered somatic variants (filtered_with_pon_vcf, ANNOVAR; VAF >= 0.04, >= 5 alt reads,
# <= 1 PoN normal), all 118 ASE-cohort tumors, no extra gnomAD filter (as Sen 2022). A and C from 62d/62e; T438/T502 exclusion as sensitivity.
# v9 Fig 3C (826/327) used all 194 tumors, PASS calls without PoN (62a); 3 ASE-cohort VEP files (T167, T172, T210) are empty.
#  A consequence counts, once per sample-gene-consequence (62b)       B missense enrichment among copy-neutral ASE cells, within-gene permutation (62f)
#  C overlap of SCNA-associated and missense-carrying PCa-ASE genes (62b)
# The v9 Fig 3D schematic and the single-gene examples are not carried over.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
P <- "/path/to/data_root/ASE"; O <- file.path(V, "pub_figs_v9"); cT <- PAL[["tumor"]]; cN <- PAL[["normal"]]; TS <- BASE / .pt
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))

cnt <- fread(file.path(V, "supp_data/S6A_consequence_counts_PoN_all.csv"))[, label := cl]
cnt[, imp := fifelse(label %like% "Missense|Inframe", "Moderate", "High")][, label := factor(label, levels = rev(label))]
ce <- fread(file.path(V, "supp_data/S6A_consequence_counts_PoN_excl_T438_T502.csv"))
note("S6A (PoN-filtered, 118 ASE-cohort tumors; once per sample-gene-consequence): %s ; excluding T438/T502: %s", cnt[, paste(sprintf("%s %d", label, N), collapse = "; ")], ce[, paste(sprintf("%s %d", cl, N), collapse = "; ")])
A <- ggplot(cnt, aes(N, label, fill = imp)) + geom_col(width = 0.65) + geom_text(aes(label = N), hjust = -0.2, size = TS, family = FONT) +
  scale_fill_manual(values = c(High = cT, Moderate = "grey55"), name = "Impact") + scale_x_continuous(expand = expansion(c(0, 0.15))) +
  labs(x = "Mutation events", y = NULL) + theme_pub() + theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))

sm <- fread(file.path(V, "supp_data/S6B_pon_permutation_summary.csv"))[set == "all" & scheme == "Within gene"]
nd <- fread(file.path(V, "supp_data/S6B_pon_permutation_null.csv"))[set == "all" & scheme == "Within gene"][, rate := 100 * mis]
obs <- 100 * sm$mis_obs
note("S6B (62f, Sen 2022 within-gene permutation, PoN-filtered, 118 tumors): copy-neutral ASE missense rate %.2f%% (%d/%d cells), expected %.2f%%, %.2f-fold, P = %.3f (%d permutations)", obs, sm$mis_n, sm$n_cells, 100 * sm$mis_null, sm$mis_fold, sm$mis_P, nrow(nd))
B <- ggplot(nd, aes(rate)) + geom_histogram(bins = 40, fill = "#FDAE6B", colour = NA) + geom_vline(xintercept = obs, colour = cT, linewidth = 0.5) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) + labs(x = "Missense rate, %", y = "Permutations") + theme_pub()

vn <- fread(file.path(V, "supp_data/S6C_venn_PoN_all.csv")); s <- setNames(vn$n, vn$set); ve <- fread(file.path(V, "supp_data/S6C_venn_PoN_excl_T438_T502.csv"))
CC <- data.table(lab = factor(c("SCNA only", "Both", "Missense only"), levels = c("SCNA only", "Both", "Missense only")), n = c(s["SCNA"] - s["overlap"], s["overlap"], s["missense"] - s["overlap"]))
note("S6C: SCNA-associated %d, missense-carrying %d, both %d (PoN-filtered, 118 tumors); excluding T438/T502: missense %d, both %d", s["SCNA"], s["missense"], s["overlap"], ve[set == "missense", n], ve[set == "overlap", n])
C <- ggplot(CC, aes(lab, n)) + geom_col(width = 0.6, fill = c("#F6C89A", "grey50", "#A9D8C5")) + geom_text(aes(label = n), vjust = -0.3, size = TS, family = FONT) +
  scale_y_continuous(expand = expansion(c(0, 0.12))) + scale_x_discrete(labels = c(`SCNA only` = "SCNA\nonly", Both = "Both", `Missense only` = "Missense\nonly")) + labs(x = NULL, y = "PCa-ASE genes") + theme_pub()

g <- (A | B | C) + plot_layout(widths = c(1.2, 1, 0.9)) + plot_annotation(tag_levels = "A") & tags
save_panel(g, file.path(O, "SuppFig_S6.pdf"), 178, 60)
writeLines(LS, file.path(O, "supp_S6_legend_stats.txt")); cat("DONE\n")
