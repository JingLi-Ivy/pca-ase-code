# 62b (2026-10-13): layout/label pass of 62 for v10 (scientific-visualization checklist); data and statistics unchanged; output pub_figs_v10.
#  S7 bottom row no longer aligned to A's long labels (blank area), 'unweighted' label; minus signs; S9B 'slope SD'; lambda -> λ.
#!/usr/bin/env Rscript
# 62: publication versions of Supplemental Figs S7 (Delta-ASE robustness), S8 (baseline/orientation sensitivity), S9 (cis-by-trans methods), S10 (cell-type ieQTL).
# Inputs are saved result files; where an estimate exists only in a log, the line is parsed with an explicit pattern and echoed
# (with its source) into supp_B_legend_stats.txt.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
O <- file.path(V, "pub_figs_v10"); cT <- PAL[["tumor"]]; cN <- PAL[["normal"]]; TS <- BASE / .pt
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
fin <- function(g, f, w, h) save_panel(g + plot_annotation(tag_levels = "A") & tags, file.path(O, f), w, h)
num <- function(l) as.numeric(regmatches(l, gregexpr("[-+]?[0-9]*\\.?[0-9]+(e-?[0-9]+)?", l))[[1]])
fp <- function(d, col = NULL) ggplot(d, aes(ratio, factor(label, levels = rev(unique(label))))) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.25, linewidth = 0.4) + labs(x = "Tumor / normal", y = NULL) + theme_pub()

l27 <- grep("^\\(A\\) filtered", readLines(file.path(V, "27_delta_ase_mapfilter.log")), value = TRUE)
l27y <- grep("weight cap", readLines(file.path(V, "27y_weight_sensitivity.log")), value = TRUE)
x <- num(sub(".*tumour/normal", "", l27)); A7 <- data.table(label = "Original weights", ratio = x[1], lo = x[2], hi = x[3], src = l27)
for (l in l27y) { v <- num(sub(".*T/N", "", l)); cap <- sub(".*weight cap +([^ ]+).*", "\\1", l); mode <- sub("^(\\w+).*", "\\1", l)
  A7 <- rbind(A7, data.table(label = sprintf("%s, cap %s", if (mode == "clean") "Overflow removed" else "Plus wrap-risk genes removed", if (cap == "none(1)") "UNW" else cap),
                             ratio = v[1], lo = v[2], hi = v[3], src = l)) }
A7[, label := sub("cap Inf", "no cap", sub("cap UNW", "unweighted", label))]
for (i in seq_len(nrow(A7))) note("S7A %s: %.3f [%.3f, %.3f]  <- %s", A7$label[i], A7$ratio[i], A7$lo[i], A7$hi[i], trimws(A7$src[i]))
A7p <- fp(A7)
ti <- fread(file.path(V, "delta_ase/tn_identity_check.csv"))
note("S7B: %d patients; tumor sites with <=1 minor-allele read: median %.3f (max %.3f); normal median %.3f", nrow(ti), median(ti$monoT), max(ti$monoT), median(ti$monoN))
B7 <- ggplot(ti, aes(100 * monoN, 100 * monoT)) + geom_abline(linetype = 2, linewidth = 0.3) + geom_point(size = 0.8, colour = PAL[["dark"]]) +
  coord_equal(xlim = c(0, 3), ylim = c(0, 3)) + labs(x = "Normal monoallelic sites, %", y = "Tumor monoallelic sites, %") + theme_pub()
rd <- function(d) rbindlist(lapply(list.files(d, full.names = TRUE), fread))[is.finite(maf_T) & is.finite(maf_N)][, .(d = mean(maf_T - maf_N)), by = pid]
mb <- merge(rd(file.path(V, "delta_mbased/one")), rd(file.path(V, "delta_mbased_phased/one")), by = "pid", suffixes = c("_pseudo", "_phased"))
note("S7C: %d patients; tumor-normal MAF difference, pseudo-phased vs phased MBASED: Spearman rho %.2f; positive in %.0f%% / %.0f%%",
     nrow(mb), cor(mb$d_pseudo, mb$d_phased, method = "s"), 100 * mean(mb$d_pseudo > 0), 100 * mean(mb$d_phased > 0))
C7 <- ggplot(mb, aes(d_pseudo, d_phased)) + geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_point(size = 0.8, colour = PAL[["dark"]]) + scale_x_continuous(breaks = c(-0.03, 0, 0.03, 0.06)) + labs(x = "T \u2212 N MAF, pseudo-phased", y = "T \u2212 N MAF, phased") + theme_pub()
tg <- function(p, t) p + labs(tag = t) + tags
save_panel(tg(A7p, "A") / wrap_elements(full = (tg(B7, "B") | tg(C7, "C"))) + plot_layout(heights = c(1.1, 1)), file.path(O, "SuppFig_S7.pdf"), 150, 125)

dl <- fread(file.path(V, "delta_ase_clean/dilution_lowbaseline_patients.csv"))
fm <- summary(lm(gain ~ dEpi, dl, weights = n))$coef
note("S8A: weak-baseline sites, %d patients; weighted regression of gain on tumor-normal epithelial difference: intercept %.4f (SE %.4f, P=%.1e), slope %.4f (P=%.1e); mean gain %.4f",
     nrow(dl), fm[1, 1], fm[1, 2], fm[1, 4], fm[2, 1], fm[2, 4], weighted.mean(dl$gain, dl$n))
A8 <- ggplot(dl, aes(dEpi, gain)) + geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") + geom_point(aes(size = n), colour = PAL[["dark"]], alpha = 0.7) +
  geom_smooth(method = "lm", mapping = aes(weight = n), colour = cT, linewidth = 0.5, se = TRUE, formula = y ~ x) +
  scale_size_area(max_size = 2, guide = "none") + labs(x = "Tumor \u2212 normal epithelial fraction", y = "Tumor \u2212 normal deviation") + theme_pub()
l40 <- grep("spec +:", readLines(file.path(V, "40_epithelial_dilution_clean.log")), value = TRUE)
B8d <- rbindlist(lapply(l40, function(l) { v <- num(sub(".*:", "", l)); data.table(base = sub(".*baseline (\\w+).*", "\\1", l), spec = sub(".*baseline \\w+ +(\\w+) spec.*", "\\1", l), ratio = v[1], lo = v[2], hi = v[3]) }))
B8d[, `:=`(base = factor(base, levels = c("low", "mid", "high"), labels = c("Weak", "Mid", "Strong")), spec = factor(spec, levels = c("low", "mid", "high"), labels = c("Low", "Mid", "High")))]
note("S8B (40_clean log): %s", B8d[, paste(sprintf("%s baseline/%s epithelial specificity %.3f [%.3f, %.3f]", base, spec, ratio, lo, hi), collapse = "; ")])
B8 <- ggplot(B8d, aes(base, ratio, colour = spec)) + geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.5), size = 0.2, linewidth = 0.4) +
  scale_colour_manual(values = c(Low = "#C6DBEF", Mid = "#6BAED6", High = "#08519C"), name = "Epithelial\nspecificity") +
  labs(x = "Normal baseline", y = "Tumor / normal") + theme_pub()
h3 <- fread(file.path(V, "delta_ase_clean/h3_trend_splits.csv")); h4 <- fread(file.path(V, "delta_ase_clean/h4_splits.csv"))[base == "all"]
lv8 <- c("Weak \u2212 strong baseline", "Tumor-preferred allele", "Normal-preferred allele")
SS <- rbind(data.table(test = lv8[1], v = h3$rdiff), data.table(test = lv8[2], v = h4[grepl("tumour", orient), ratio - 1]),
            data.table(test = lv8[3], v = h4[grepl("normal", orient), ratio - 1]))
note("S8C: across 20 random splits: %s", SS[, .(r = sprintf("%.3f-%.3f", min(v), max(v))), by = test][, paste(test, r, collapse = "; ")])
C8 <- ggplot(SS, aes(v, factor(test, levels = rev(lv8)))) + geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) + geom_jitter(height = 0.15, width = 0, size = 0.7, colour = PAL[["dark"]]) +
  labs(x = "Estimate across 20 splits", y = NULL) + theme_pub()
fin((A8 | B8 | C8) + plot_layout(widths = c(1, 1.1, 1.1)), "SuppFig_S8.pdf", 178, 58)

sim <- fread(file.path(V, "reliability_bias/simulation_reps.csv"))
SM <- melt(sim, id.vars = c("r", "tissue", "rep"), measure.vars = c("fz", "slz"))[, .(m = mean(value), s = sd(value)), by = .(r, variable)]
SM[, stat := factor(variable, levels = c("fz", "slz"), labels = c("Correlation difference", "Slope difference"))]
note("S9A (49): simulated statistic with no cis-by-trans (equal true slopes), by true |a|-TF correlation: %s", SM[, paste(sprintf("%s r=%.2f: %.3f", stat, r, m), collapse = "; ")])
A9 <- ggplot(SM, aes(r, m, colour = stat)) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) + geom_line(linewidth = 0.4) +
  geom_pointrange(aes(ymin = m - s, ymax = m + s), size = 0.2, linewidth = 0.3) +
  scale_colour_manual(values = c(`Correlation difference` = PAL[["accent"]], `Slope difference` = cN), name = NULL) +
  labs(x = "True |a|-TF correlation", y = "Spurious interaction, z") + theme_pub() + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
l50 <- readLines(file.path(V, "50_h5_bias_expectation.log")); i <- grep("expected_H5 +observed_H5_here", l50)
B9d <- rbindlist(lapply(l50[(i + 2):(i + 4)], function(l) { v <- num(l); data.table(sb = v[2], expected = v[5], observed = v[6]) }))
note("S9B (50 log): expected H5 under no cis-by-trans %s; observed on the same 30-TF subsample %.4f", paste(sprintf("%.4f (sb %.2f)", B9d$expected, B9d$sb), collapse = ", "), B9d$observed[1])
B9 <- ggplot(rbind(B9d[, .(lab = sprintf("Expected, slope SD %.2f", sb), v = expected)], data.table(lab = "Observed", v = B9d$observed[1])), aes(v, factor(lab, levels = rev(lab)))) +
  geom_col(width = 0.6, fill = c(rep("grey65", 3), cT)) + labs(x = "Tumor \u2212 normal coupling, z", y = NULL) + theme_pub()
cm <- fread(file.path(V, "cistrans_slope_26e/tn_slope_by_covariate_model.csv"))
cm[, label := c(`ai ~ ld` = "Depth", `ai ~ ld + epi` = "+ epithelial score", `ai ~ ld + pro` = "+ proliferation", `ai ~ ld + epi + pro` = "+ both", `ai ~ ld + bp_epi` = "+ BayesPrism epithelial")[model]]
note("S9C (26e): %s", cm[, paste(sprintf("%s %.4f [%.4f, %.4f]", label, delta, delta_lo, delta_hi), collapse = "; ")])
C9 <- ggplot(cm, aes(delta, factor(label, levels = rev(label)))) + geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_pointrange(aes(xmin = delta_lo, xmax = delta_hi), size = 0.25, linewidth = 0.4) + labs(x = "Tumor \u2212 normal coupling, z", y = NULL) + theme_pub()
l28 <- grep("M1b lambda", readLines(file.path(V, "28f_ieqtl_specificity_gtex_paired.log")), value = TRUE)
D9d <- rbindlist(lapply(l28, function(l) { rg <- as.numeric(strsplit(sub(".*range ([0-9.]+-[0-9.]+);.*", "\\1", l), "-")[[1]])
  data.table(tissue = ifelse(grepl("^\\[T\\]", l), "Tumor", "Normal"), tf = as.numeric(sub(".*TFs ([0-9.]+) .*", "\\1", l)),
  rnd = as.numeric(sub(".*non-TF genes ([0-9.]+) .*", "\\1", l)), lo = rg[1], hi = rg[2], src = l) }))
for (k in seq_len(nrow(D9d))) note("S9D %s: TF lambda %.3f; matched random non-TF genes %.3f (range %.3f-%.3f, 1000 draws)  <- %s", D9d$tissue[k], D9d$tf[k], D9d$rnd[k], D9d$lo[k], D9d$hi[k], trimws(D9d$src[k]))
D9 <- ggplot(D9d, aes(y = tissue)) + geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, colour = "grey55", orientation = "y", linewidth = 0.5) +
  geom_point(aes(x = rnd), colour = "grey55", size = 1.2) + geom_point(aes(x = tf), colour = cT, size = 1.6, shape = 18) +
  labs(x = "Interaction eQTL inflation, \u03bb", y = NULL) + theme_pub()
fin((A9 | B9) / (C9 | D9), "SuppFig_S9.pdf", 178, 105)

ct <- fread(file.path(V, "celltype_ieqtl/celltype_ieqtl_summary.csv")); pr <- fread(file.path(V, "celltype_ieqtl/celltype_ieqtl_primary_1000perm.csv"))
ct <- ct[leads == "GTEx" & purity == FALSE]
ct <- merge(ct[, !c("null_lambda", "null_lambda_max", "P_lambda")], rbind(pr[, .(tissue, cell, null_lambda, null_lambda_max, P_lambda)],
            ct[cell != "Epithelial_cells", .(tissue, cell, null_lambda, null_lambda_max, P_lambda)]), by = c("tissue", "cell"))
ct[, `:=`(tis = factor(tissue, levels = c("N", "T"), labels = c("Normal", "Tumor")), cel = factor(cell, levels = c("Epithelial_cells", "Fibroblast", "SMC"), labels = c("Epithelial", "Fibroblast", "Smooth muscle")))]
note("S10: %s", ct[order(tis, cel), paste(sprintf("%s %s lambda %.2f (null max %.2f; P=%.3f; %d genes FDR<0.05)", tis, cel, lambda, null_lambda_max, P_lambda, fdr05), collapse = "; ")])
note("S10: epithelial P from 1,000 permutations; fibroblast and smooth muscle from 100 permutations")
A10 <- ggplot(ct, aes(lambda, cel, colour = tis)) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_linerange(aes(xmin = null_lambda, xmax = null_lambda_max, group = tis), position = position_dodge(0.5), colour = "grey80", linewidth = 1.5) +
  geom_point(position = position_dodge(0.5), size = 1.6) + scale_colour_manual(values = c(Normal = cN, Tumor = cT), name = NULL) +
  labs(x = "Interaction inflation, \u03bb", y = NULL) + theme_pub() + theme(legend.position = "top")
B10 <- ggplot(ct, aes(signed_mean_t, cel, colour = tis)) + geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_pointrange(aes(xmin = signed_lo, xmax = signed_hi), position = position_dodge(0.5), size = 0.25, linewidth = 0.4) +
  scale_colour_manual(values = c(Normal = cN, Tumor = cT), guide = "none") + labs(x = "Signed interaction t", y = NULL) + theme_pub()
fin(A10 | B10, "SuppFig_S10.pdf", 140, 55)
writeLines(LS, file.path(O, "supp_B_legend_stats.txt")); cat("DONE\n")
