#!/usr/bin/env Rscript
# 61: publication versions of Supplemental Figs S1 (ASE calling QC), S2 (transcriptome-wide T vs N), S3 (GO),
#     S4 (QTL enrichment, eQTL direction), S5 (methylation triangulation), S15 (global imbalance vs Gleason).
# Inputs are saved results of the original analyses (paths below); S1 A-C plot data recomputed by 61a with 07's seed and procedure.
# Style: fig_theme.R (Arial 8 pt, final size). Legend numbers -> supp_A_legend_stats.txt.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
P <- "/path/to/data_root/ASE"; SD <- file.path(V, "supp_data"); O <- file.path(V, "pub_figs_v9"); dir.create(O, showWarnings = FALSE)
cT <- PAL[["tumor"]]; cN <- PAL[["normal"]]; TS <- BASE / .pt
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
fin <- function(g, f, w, h, tag = TRUE) save_panel(if (tag) g + plot_annotation(tag_levels = "A") & tags else g, file.path(O, f), w, h)

s1 <- fread(file.path(SD, "S1_summary.csv")); h1 <- fread(file.path(SD, "S1A_depth_hist.csv")); h2 <- fread(file.path(SD, "S1B_refratio_hist.csv"))
qs <- fread(file.path(SD, "S1C_summary.csv")); qp <- fread(file.path(SD, "S1C_null_p.csv")); ds <- fread(file.path(P, "code/ase_qc_results/ASE_depth_sensitivity.csv"))
note("S1: %s heterozygous SNP observations in %d tumors; median depth %g; background sites (depth>=10, FDR>=0.1) %s, ref fraction median %.3f mean %.4f",
     format(s1$n_obs, big.mark = ","), s1$n_samples, s1$median_depth, format(s1$bg_n, big.mark = ","), s1$bg_median, s1$bg_mean)
note("S1C: %d background sites; overdispersion rho %.4f; lambda beta-binomial %.2f, binomial %.2f", qs$n, qs$rho, qs$lambda_betabinom, qs$lambda_binom)
note("S1D: significant ASE genes by minimum depth: %s", ds[, paste(sprintf(">=%d: %d", min_depth, n_gene_sig), collapse = "; ")])
A <- ggplot(h1[depth <= 200], aes(depth, N / 1e6)) + geom_col(fill = "#4393C3", width = 1) + geom_vline(xintercept = 10, linetype = 2, linewidth = 0.3) +
  scale_x_continuous(limits = c(0, 100)) + labs(x = "Allelic read depth", y = "SNP observations, million") + theme_pub()
h2b <- h2[, .(N = sum(N)), by = .(bin = 2 * (bin %/% 2))]
B <- ggplot(h2b, aes((bin + 1) / 100, N / 1e3)) + geom_col(fill = "grey55", width = 0.02) + geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3) +
  labs(x = "Reference-allele fraction", y = "Background sites, thousand") + theme_pub()
qq <- rbind(data.table(model = "Beta-binomial", p = sort(qp$p_betabinom)), data.table(model = "Binomial", p = sort(qp$p_binom)))
qq[, exp := -log10(ppoints(.N)), by = model][, obs := -log10(pmax(p, 1e-300))]
C <- ggplot(qq, aes(exp, obs, colour = model)) + geom_abline(linetype = 2, linewidth = 0.3) + geom_point(size = 0.3) +
  scale_colour_manual(values = c(`Beta-binomial` = cN, Binomial = cT), name = NULL) + labs(x = "Expected -log10 P", y = "Observed -log10 P") +
  theme_pub() + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1)) + guides(colour = guide_legend(override.aes = list(size = 1.5)))
D <- ggplot(ds, aes(factor(min_depth), n_gene_sig / 1e3, group = 1)) + geom_line(colour = cN, linewidth = 0.4) + geom_point(colour = cN, size = 1.2) +
  labs(x = "Minimum read depth", y = "Genes with significant ASE, thousand") + theme_pub()
fin(A | B | C | D, "SuppFig_S1.pdf", 178, 50)

S <- fread(file.path(P, "figure_freeze/Fig1/data/Fig1D_threeset_burden_data.csv"))[set == "All genes"]
w <- wilcox.test(S$Tumor, S$Normal, paired = TRUE)
note("S2A: %d genes with data in both tissues; median |a| tumor %.3f vs normal %.3f; %.1f%% genes higher in tumor; paired Wilcoxon P=%.2g",
     nrow(S), median(S$Tumor), median(S$Normal), 100 * mean(S$Tumor > S$Normal), w$p.value)
A2 <- ggplot(S, aes(Normal, Tumor)) + geom_bin2d(bins = 60) + geom_abline(linetype = 2, linewidth = 0.3) +
  scale_fill_gradient(low = "grey90", high = "black", trans = "log10", name = "Genes") + coord_equal() +
  labs(x = "Normal mean |a|", y = "Tumor mean |a|") + theme_pub()
pca <- fread(file.path(P, "code/PCa_ASE_genes.csv"))
fp <- fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls")
tc <- grep("^T\\d+_WTS$", names(fp), value = TRUE); fp[, medFPKM := apply(.SD, 1, median, na.rm = TRUE), .SDcols = tc]
m <- merge(pca[, .(gene_symbol, pct = as.numeric(pct_sig_Tumor))], fp[, .(gene_symbol = gene_name, medFPKM)], by = "gene_symbol")[is.finite(medFPKM) & medFPKM > 0 & is.finite(pct)]
m[, grp := cut(medFPKM, quantile(medFPKM, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("Low", "Medium", "High"))]
kw <- kruskal.test(pct ~ grp, m)
note("S2B: %d PCa-ASE genes with expression; median %% tumors with ASE by tertile %s; Kruskal-Wallis chi2=%.1f, df=2, P=%.2g",
     nrow(m), m[, .(v = sprintf("%.1f", median(pct))), by = grp][order(grp), paste(grp, v, collapse = "/")], kw$statistic, kw$p.value)
B2 <- ggplot(m, aes(grp, pct)) + geom_boxplot(width = 0.6, linewidth = 0.3, outlier.size = 0.2, fill = "grey85") +
  coord_cartesian(ylim = c(0, quantile(m$pct, 0.99))) + labs(x = "Tumor expression tertile", y = "Tumors with significant ASE, %") + theme_pub()
fin((A2 | B2) + plot_layout(widths = c(1.2, 1)), "SuppFig_S2.pdf", 120, 55)

go <- fread(file.path(P, "code/enrichment_results/GO_enrichment_BP_significant.csv"))[order(Fisher_pvalue)][1:20]
note("S3: top 20 of %d GO BP terms with FDR<0.1 (topGO Fisher; background 28,851 genes with heterozygous-SNP coverage); FDR range %.1e-%.1e",
     nrow(fread(file.path(P, "code/enrichment_results/GO_enrichment_BP_significant.csv"))), min(go$FDR), max(go$FDR))
go[, Term := factor(Term, levels = rev(Term))]
G3 <- ggplot(go, aes(Enrichment_ratio, Term, size = Significant_num, colour = -log10(FDR))) + geom_point() +
  scale_colour_gradient(low = "#92C5DE", high = cT, name = "-log10 FDR") + scale_size_area(max_size = 3, name = "Genes") +
  labs(x = "Enrichment ratio", y = NULL) + theme_pub()
fin(G3, "SuppFig_S3.pdf", 150, 95, tag = FALSE)

d6 <- fread(file.path(P, "figures/pub_figs/S6D_data.csv")); f6 <- fread(file.path(P, "figures/pub_figs/S6F_data.csv"))
for (i in seq_len(nrow(d6))) note("S4A %s vs %s: OR %.2f [%.2f, %.2f], P=%.2g", d6$qtl[i], d6$bg[i], d6$OR[i], d6$lo[i], d6$hi[i], d6$p[i])
note("S4B: OR by expression tertile x MAF bin: %s", f6[, paste(sprintf("%s/%s %.2f", etert, mafb, OR), collapse = "; ")])
A4 <- ggplot(d6, aes(OR, bg, colour = qtl)) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.25, linewidth = 0.4, position = position_dodge(0.5)) +
  scale_colour_manual(values = c(eQTL = cN, sQTL = PAL[["accent"]]), name = NULL) + labs(x = "Odds ratio", y = NULL) + theme_pub() + theme(legend.position = "top")
f6[, etert := factor(etert, levels = c("Low", "Mid", "High"))]
B4 <- ggplot(f6, aes(mafb, OR, colour = etert, group = etert)) + geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3) + geom_line(linewidth = 0.4) + geom_point(size = 1) +
  scale_colour_manual(values = c(Low = "#9ECAE1", Mid = "#4292C6", High = "#08519C"), name = "Expression") + labs(x = "Minor-allele frequency", y = "eQTL odds ratio") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
cm <- fread(file.path(P, "figures/eqtl_concordance_full/eqtl_full_matched_obs.csv"))
sn <- cm[, .(ref_n = sum(ref_n), alt_n = sum(alt_n), eqtl_high = eqtl_high[1], ref = ref[1], alt = alt[1]), by = .(chrom, pos)][ref_n != alt_n]
sn[, conc := eqtl_high == fifelse(ref_n > alt_n, ref, alt)]
cc <- data.table(level = c("Per\nobservation", "Per\nSNP"), k = c(sum(cm$concordant), sum(sn$conc)), n = c(nrow(cm), nrow(sn)))
cc[, `:=`(f = k / n, P = mapply(function(k, n) binom.test(k, n)$p.value, k, n))]
note("S4C: eQTL-ASE direction concordance: %s (%d genes)", cc[, paste(sprintf("%s %d/%d = %.1f%%, binomial P=%.1e", level, k, n, 100 * f, P), collapse = "; ")], uniqueN(cm$gene))
C4 <- ggplot(cc, aes(level, f)) + geom_col(width = 0.55, fill = "grey45") + geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * f)), vjust = -0.3, size = TS, family = FONT) + scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
  labs(x = NULL, y = "Direction concordance") + theme_pub()
fin((A4 | B4 | C4) + plot_layout(widths = c(1.3, 1.2, 0.7)), "SuppFig_S4.pdf", 178, 58)

mt <- fread(file.path(P, "code/methylation_triangulation_results/methylation_triangulation_per_gene.csv"))
note("S5: %d genes with promoter methylation; methylation-ASE FDR<0.1 %d; triangulated %d", nrow(mt), mt[fdr_meth_ase < 0.1, .N], mt[triangulated %in% TRUE, .N])
mt[, cl := fifelse(triangulated %in% TRUE, "Coherent", "Other")]
lab <- mt[triangulated %in% TRUE][order(fdr_meth_ase)][1:10]
S5 <- ggplot(mt[order(cl == "Coherent")], aes(rho_meth_expr, rho_meth_ase, colour = cl)) + geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") + geom_point(size = 0.6) +
  geom_text(data = lab, aes(label = gene_symbol), colour = "black", size = TS, family = FONT, fontface = "italic", vjust = -0.6, check_overlap = TRUE) +
  scale_colour_manual(values = c(Coherent = cT, Other = "grey75"), name = NULL) +
  labs(x = "Methylation-expression correlation", y = "Methylation-ASE correlation") + theme_pub() + theme(legend.position = "top")
fin(S5 + scale_y_continuous(expand = expansion(c(0.05, 0.12))) + coord_cartesian(clip = "off"), "SuppFig_S5.pdf", 85, 80, tag = FALSE)

rob <- fread(file.path(P, "figures/ase_clinical/ase_clinical_robust_data.csv"))[, sid := as.character(sid)]
tcn <- fread(file.path(P, "figures/ase_clinical/ase_gleason_true_copyneutral.csv"))[, sid := as.character(sid)]
dA <- rob[is.finite(gleason) & is.finite(imb_all)]; sp <- cor.test(dA$imb_all, dA$gleason, method = "spearman", exact = FALSE)
d <- merge(rob[, .(sid, gleason, imb_all)], tcn[, .(sid, imb_baldiploid)], by = "sid")[is.finite(gleason) & is.finite(imb_all) & is.finite(imb_baldiploid)]
dm <- melt(d[, .(sid, grade = factor(fifelse(gleason >= 8, "8-10", "6-7"), levels = c("6-7", "8-10")), All = imb_all, `Copy-neutral` = imb_baldiploid)], id.vars = c("sid", "grade"))
pv <- dm[, .(P = wilcox.test(value ~ grade)$p.value, n_lo = sum(grade == "6-7"), n_hi = sum(grade == "8-10")), by = variable]
note("S15A: n=%d tumors; Spearman rho=%.2f, P=%.1e", nrow(dA), sp$estimate, sp$p.value)
note("S15B: %s", pv[, paste(sprintf("%s sites: Gleason 6-7 n=%d vs 8-10 n=%d, Mann-Whitney P=%.1e", variable, n_lo, n_hi, P), collapse = "; ")])
A15 <- ggplot(dA, aes(factor(gleason), imb_all)) + geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3, fill = "grey85") +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.5) + labs(x = "Gleason score", y = "Mean allelic imbalance") + theme_pub()
B15 <- ggplot(dm, aes(variable, value, fill = grade)) + geom_boxplot(outlier.size = 0.2, width = 0.65, linewidth = 0.3) +
  scale_fill_manual(values = c(`6-7` = cN, `8-10` = cT), name = "Gleason") + labs(x = "Sites", y = "Mean allelic imbalance") + theme_pub() + theme(legend.position = "top")
fin(A15 | B15, "SuppFig_S15.pdf", 120, 55)
writeLines(LS, file.path(O, "supp_A_legend_stats.txt")); cat("DONE\n")
