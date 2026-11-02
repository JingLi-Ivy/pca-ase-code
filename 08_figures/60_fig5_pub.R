#!/usr/bin/env Rscript
# 60: publication version of main Fig 5 — testing transcription-factor dosage amplification at two scales (R5, R6).
#  A H5: tumor-normal cis-by-trans slope difference vs 1,000 whole-patient label permutations (26e, model ai ~ ld)
#  B prostate lineage / driver TF sets vs expression-matched TF sets (52, primary statistic = tumor-normal difference)
#  C KIF18B-KIF18B-DT bidirectional promoter: exons, rs9914422, reporter fragment, ENCODE cCREs (locus_fig inputs of 22)
#  D ADASTRA v6.1 allele-specific binding at rs9914422 (FDR < 0.05, plus FOXM1)
#  E phased, DNA-corrected KIF18B allelic output vs FOXM1 (kif18b_phase/phased_patient_level.csv, as 25)
#  F dual-luciferase reporter, promoter (pGL4.10) and enhancer (pGL4.23) constructs, empty vector vs FOXM1
#  G siFOXM1 qRT-PCR in PC-3 (FOXM1, KIF18B, NEK2, TOP2A)
#  (H revised model: schematic, finished in Illustrator)
# F/G read wetlab/ASE_raw_data.xlsx (copy of submission/ASE-Raw -data.xlsx). Block assignment of the luciferase sheet:
# columns B-D/G-I = promoter (allele main effect P = 0.024 as in v9), M-O/R-T = enhancer (P = 3.3e-4).
# NOTE 2026-10-07: enhancer replicate values repeat (N3=N5 2.35481, S4=S5 3.478192, T4=T5 2.88166). Checked against the
#   plate-reader records on 2026-10-28: the two wells gave identical readings, so the ratios are identical (not duplicated entries).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(jsonlite); library(readxl)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
O <- file.path(V, "pub_figs_v9"); dir.create(O, showWarnings = FALSE)
cT <- PAL[["tumor"]]; cN <- PAL[["normal"]]
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
TS <- BASE / .pt

TN <- fread(file.path(V, "cistrans_slope_26e/tn_slope_by_covariate_model.csv"))[model == "ai ~ ld"]
NL <- fread(file.path(V, "cistrans_slope_26e/tn_slope_strict_null.csv"))[model == "r0"]
pp <- (1 + sum(abs(NL$delta - mean(NL$delta)) >= abs(TN$delta - mean(NL$delta)))) / (1 + nrow(NL))
note("Fig5A: observed %.4f [%.4f, %.4f]; %d permutations, null mean %.4f SD %.4f; %.1f SD; two-sided P=%.3f",
     TN$delta, TN$delta_lo, TN$delta_hi, nrow(NL), mean(NL$delta), sd(NL$delta), (TN$delta - mean(NL$delta)) / sd(NL$delta), pp)
A5 <- ggplot(NL, aes(delta)) + geom_histogram(bins = 30, fill = "grey75", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = TN$delta, colour = cT, linewidth = 0.6) + geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) + labs(x = "Tumor - normal coupling, z", y = "Permutations") + theme_pub()

TFS <- fread(file.path(V, "prostate_tf_set/set_vs_matched_tf.csv"))[model == "r0" & stat == "delta"]
TFS[, lab := factor(fifelse(set == "LIN", "Lineage TFs", "Driver TFs"), levels = c("Driver TFs", "Lineage TFs"))]
for (i in seq_len(nrow(TFS))) note("Fig5B %s (%s): observed %.4f; matched null %.4f +- %.4f; P=%.2f; all-TF mean %.4f", TFS$lab[i], TFS$members[i], TFS$obs[i], TFS$matched_null_mean[i], TFS$matched_null_sd[i], TFS$P[i], TFS$all_TF_mean[i])
B5 <- ggplot(TFS, aes(y = lab)) + geom_vline(xintercept = TFS$all_TF_mean[1], linetype = 2, linewidth = 0.3) +
  geom_errorbar(aes(xmin = matched_null_mean - 1.96 * matched_null_sd, xmax = matched_null_mean + 1.96 * matched_null_sd), width = 0.25, orientation = "y", colour = "grey60", linewidth = 0.5) +
  geom_point(aes(x = obs), colour = cT, size = 1.8) + scale_x_continuous(breaks = c(-0.02, 0.02, 0.06)) + labs(x = "Tumor - normal coupling, z", y = NULL) + theme_pub()

L <- file.path(V, "locus_fig"); RS <- 44947872; FR <- c(44947189, 44948553); ZM <- c(44946900, 44949000)
ex <- fread(file.path(L, "gene_exons.tsv"), header = FALSE, col.names = c("gene", "tx", "type", "start", "end", "strand"))[type == "exon" & end >= ZM[1] & start <= ZM[2]]
ex[, `:=`(s = pmax(start, ZM[1]), e = pmin(end, ZM[2]))]; ex <- unique(ex[, .(gene, s, e)])
yg <- c(KIF18B = 3.3, `KIF18B-DT` = 2.2)
cc <- as.data.table(fromJSON(file.path(L, "ccre.json"))); cc[, cls := get(if ("ucscLabel" %in% names(cc)) "ucscLabel" else "encodeLabel")]
cc <- cc[chromEnd >= ZM[1] & chromStart <= ZM[2]][, `:=`(s = pmax(chromStart, ZM[1]), e = pmin(chromEnd, ZM[2]))]
cc[, cl2 := fifelse(cls %in% c("prom", "PLS"), "Promoter-like", fifelse(cls %in% c("enhP", "pELS", "enhD", "dELS"), "Enhancer-like", "Other"))]
note("Fig5C: window chr17:%s-%s (hg38); rs9914422 at %s; reporter fragment %s-%s (%d bp); cCREs %s",
     format(ZM[1], big.mark = ","), format(ZM[2], big.mark = ","), format(RS, big.mark = ","), format(FR[1], big.mark = ","), format(FR[2], big.mark = ","), diff(FR) + 1, paste(unique(cc$cls), collapse = ","))
C5 <- ggplot() +
  geom_segment(data = data.table(gene = names(yg), y = yg), aes(x = ZM[1], xend = ZM[2], y = y, yend = y), colour = "grey70", linewidth = 0.3) +
  geom_rect(data = ex, aes(xmin = s, xmax = e, ymin = yg[gene] - 0.25, ymax = yg[gene] + 0.25, fill = gene)) +
  annotate("text", x = ZM[1] + 50, y = c(3.8, 2.7), label = c("KIF18B, minus strand", "KIF18B-DT, plus strand"), hjust = 0, size = TS, family = FONT, fontface = "italic") +
  annotate("rect", xmin = FR[1], xmax = FR[2], ymin = 1.15, ymax = 1.4, fill = "grey55") +
  annotate("text", x = ZM[1] + 50, y = 1.65, label = "Reporter fragment", hjust = 0, size = TS, family = FONT) +
  geom_rect(data = cc, aes(xmin = s, xmax = e, ymin = 0.2, ymax = 0.6, fill = cl2)) +
  annotate("text", x = ZM[1] + 50, y = 0.82, label = "ENCODE cCRE", hjust = 0, size = TS, family = FONT) +
  geom_vline(xintercept = RS, colour = cT, linewidth = 0.4) + annotate("text", x = RS + 40, y = 4.3, label = "rs9914422", hjust = 0, size = TS, family = FONT, colour = cT) +
  scale_fill_manual(values = c(KIF18B = cN, `KIF18B-DT` = "#5AAE61", `Promoter-like` = "#E41A1C", `Enhancer-like` = "#FFA700", Other = "grey70"), guide = "none") +
  scale_x_continuous(limits = ZM, expand = c(0, 0), breaks = c(44947000, 44948000), labels = c("44.947", "44.948")) +
  scale_y_continuous(limits = c(0, 4.5)) + coord_cartesian(clip = "off") + labs(x = "chr17, Mb", y = NULL) + theme_pub() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

ad <- fread(file.path(L, "adastra_rs9914422.tsv"))[minfdr < 0.05 | tf == "FOXM1"][order(minfdr)]
ad[, tf := factor(tf, levels = rev(tf))]
note("Fig5D: ADASTRA FDR<0.05 factors: %s ; FOXM1 FDR %.3f", ad[minfdr < 0.05, paste(sprintf("%s(%s) %.1e", tf, pref, minfdr), collapse = ", ")], ad[tf == "FOXM1", minfdr])
D5 <- ggplot(ad, aes(-log10(minfdr), tf, fill = pref)) + geom_col(width = 0.65) + geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = c(G = "grey45", C = cT), name = "Preferred allele") + scale_x_continuous(expand = expansion(c(0, 0.05))) +
  labs(x = "-log10 FDR, allele-specific binding", y = NULL) + theme_pub() + theme(legend.position = "top", legend.margin = margin(0, 0, -2, 0, "mm"))

Ph <- fread(file.path(V, "kif18b_phase/phased_patient_level.csv"))[rnaC + rnaG > 0]
Ph[, n := rnaC + rnaG]; Ph[, c("lo", "hi") := as.list(binom.test(rnaC, n)$conf.int), by = pid]
Ph[, state := fifelse(tcn == 2 & lcn == 1, "Balanced", "Copy-number altered")]
ct <- cor.test(Ph$A_reg, Ph$FOXM1, method = "s", exact = FALSE)
note("Fig5E: %d tumors, %d phased reads; pooled C fraction %.3f; DNA-corrected A_reg vs FOXM1 Spearman rho %.2f, P=%.2f (beta-binomial: C fraction 0.477 [0.39, 0.58], P=0.77; FOXM1 slope P=0.60, balanced P=0.89, from 08b)",
     nrow(Ph), sum(Ph$n), sum(Ph$rnaC) / sum(Ph$n), ct$estimate, ct$p.value)
E5 <- ggplot(Ph, aes(FOXM1, rna_frac_C)) + annotate("rect", xmin = min(Ph$FOXM1) * 0.9, xmax = max(Ph$FOXM1) * 1.1, ymin = 0.39, ymax = 0.58, fill = "grey88") +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) + geom_errorbar(aes(ymin = lo, ymax = hi, colour = state), width = 0, alpha = 0.4, linewidth = 0.3) +
  geom_point(aes(size = n, colour = state)) + scale_x_log10() + scale_size_area(max_size = 2.5, name = "Reads") +
  scale_colour_manual(values = c(Balanced = "black", `Copy-number altered` = PAL[["accent"]]), name = NULL) +
  labs(x = "FOXM1, FPKM", y = "C-haplotype fraction") + theme_pub() + theme(legend.position = "right")

XL <- file.path(V, "wetlab/ASE_raw_data.xlsx")
lu <- as.data.frame(read_excel(XL, sheet = "Dual luciferase", col_names = FALSE, .name_repair = "minimal"))
val <- function(r, c) as.numeric(lu[r, c])
blk <- function(cols, construct, fox) rbindlist(lapply(seq_along(cols), function(j) data.table(construct = construct, fox = fox, allele = c("Control", "C", "G")[j], y = sapply(2:4, function(r) val(r, cols[j])))))
LU <- rbind(blk(2:4, "Promoter", "Vector"), blk(7:9, "Promoter", "FOXM1"), blk(13:15, "Enhancer", "Vector"), blk(18:20, "Enhancer", "FOXM1"))
LU[, `:=`(allele = factor(allele, levels = c("Control", "C", "G")), fox = factor(fox, levels = c("Vector", "FOXM1")), construct = factor(construct, levels = c("Promoter", "Enhancer")))]
LM <- LU[, .(m = mean(y), s = sd(y)), by = .(construct, fox, allele)]
for (cn in c("Promoter", "Enhancer")) { d <- LU[construct == cn & allele != "Control"]; a <- anova(lm(log(y) ~ allele * fox, d))
  r <- LM[construct == cn]; note("Fig5F %s: C/G vector %.2f, C/G FOXM1 %.2f; two-way ANOVA on log: allele P=%.2g, FOXM1 P=%.2g, interaction P=%.2g (n=3 technical replicates, single experiment)",
    cn, r[fox == "Vector" & allele == "C", m] / r[fox == "Vector" & allele == "G", m], r[fox == "FOXM1" & allele == "C", m] / r[fox == "FOXM1" & allele == "G", m], a[1, 5], a[2, 5], a[3, 5]) }
F5 <- ggplot(LM, aes(allele, m, fill = allele)) + geom_col(width = 0.65) + geom_errorbar(aes(ymin = m - s, ymax = m + s), width = 0.2, linewidth = 0.3) +
  geom_point(data = LU, aes(allele, y), size = 0.5, position = position_jitter(width = 0.1, seed = 60)) +
  facet_grid(. ~ construct + fox) + scale_fill_manual(values = c(Control = "grey85", C = cT, G = "grey45"), guide = "none") + scale_x_discrete(labels = c(Control = "E", C = "C", G = "G")) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) + labs(x = NULL, y = "Relative luciferase activity") + theme_pub() +
  theme(panel.spacing = unit(1, "mm"))

qp <- as.data.frame(read_excel(XL, sheet = "qPCR2", range = "A1:K5", col_names = FALSE, .name_repair = "minimal"))
QP <- rbindlist(lapply(2:5, function(r) data.table(gene = as.character(qp[r, 2]), cond = rep(c("siNC", "siFOXM1-1", "siFOXM1-2"), each = 3), y = as.numeric(unlist(qp[r, 3:11])))))
QP[, `:=`(gene = factor(gene, levels = c("FOXM1", "KIF18B", "NEK2", "TOP2A")), cond = factor(cond, levels = c("siNC", "siFOXM1-1", "siFOXM1-2")))]
QM <- QP[, .(m = mean(y), s = sd(y)), by = .(gene, cond)]
note("Fig5G: mean relative expression (siFOXM1-1 / siFOXM1-2): %s", QM[cond != "siNC", .(v = paste(sprintf("%.2f", m), collapse = "/")), by = gene][, paste(gene, v, collapse = "; ")])
G5 <- ggplot(QM, aes(gene, m, fill = cond)) + geom_col(width = 0.7, position = position_dodge(0.75)) +
  geom_errorbar(aes(ymin = m - s, ymax = m + s), width = 0.2, linewidth = 0.3, position = position_dodge(0.75)) +
  scale_fill_manual(values = c(siNC = "grey80", `siFOXM1-1` = "#4393C3", `siFOXM1-2` = "#2166AC"), name = NULL) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) + labs(x = NULL, y = "Relative mRNA") + theme_pub() +
  theme(axis.text.x = element_text(face = "italic"), legend.position = "top", legend.margin = margin(0, 0, -2, 0, "mm"))

g <- ((A5 | B5 | C5) + plot_layout(widths = c(1, 1, 1.3))) / ((D5 | E5) + plot_layout(widths = c(1, 1.4))) / ((F5 | G5) + plot_layout(widths = c(2, 1))) +
  plot_annotation(tag_levels = "A") & tags
save_panel(g, file.path(O, "Fig5_AG.pdf"), 178, 170)
writeLines(LS, file.path(O, "fig5_legend_stats.txt")); cat("DONE\n")
