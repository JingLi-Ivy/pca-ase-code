#!/usr/bin/env Rscript
# 57: publication versions of main Fig 3 (R3, H1) and Fig 4 (R4, H2-H4), redrawn from the CLEAN results
# (56: overflow rows and wrap-risk genes removed; weights capped at 100 reads). Replaces 47 for the manuscript.
# Every number is read from a result file (no hand-typed estimates); legend statistics are written to fig3_fig4_legend_stats.txt.
# Style (fig_theme.R): final size 178 mm wide, Arial 8 pt, no in-panel titles, American spelling.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)}); set.seed(57)
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
D <- file.path(V, "delta_ase_clean"); E <- file.path(V, "delta_eqtl_clean"); O <- file.path(V, "pub_figs_v9"); dir.create(O, showWarnings = FALSE)
cT <- PAL[["tumor"]]; cN <- PAL[["normal"]]; TXT <- BASE / .pt
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))

P <- fread(file.path(D, "mapfilter_A_per_patient.csv"))
bt <- replicate(2000, { p <- P[sample(.N, replace = TRUE)]; p[, weighted.mean(sT, n)] / p[, weighted.mean(sN, n)] })
note("Fig3A: n=%d patients; tumor/normal %.3f [%.3f, %.3f]; %.0f%% sT>sN; paired Wilcoxon P=%.1g", nrow(P),
     P[, weighted.mean(sT, n)] / P[, weighted.mean(sN, n)], quantile(bt, .025), quantile(bt, .975), 100 * P[, mean(sT > sN)], wilcox.test(P$sT, P$sN, paired = TRUE)$p.value)
lim <- range(c(P$sN, P$sT)); lim <- lim + c(-1, 1) * 0.05 * diff(lim)
A3 <- ggplot(P, aes(sN, sT)) + geom_abline(linetype = 2, colour = "grey50", linewidth = 0.3) + geom_point(size = 0.8, colour = PAL[["dark"]], alpha = 0.8) +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) + labs(x = "Normal allelic deviation", y = "Tumor allelic deviation") + theme_pub()

B <- fread(file.path(D, "mapfilter_B_per_snp.csv.gz"))[is.finite(p)]
B[, cl := factor(fifelse(sig %in% TRUE, fifelse(b > 0, "Tumor stronger", "Normal stronger"), "NS"), levels = c("Tumor stronger", "Normal stronger", "NS"))]
note("Fig3B: %s SNPs tested; FDR<0.05 & |diff|>=0.05: %d tumor-stronger, %d normal-stronger", format(nrow(B), big.mark = ","), B[cl == "Tumor stronger", .N], B[cl == "Normal stronger", .N])
B3 <- ggplot(B[order(cl != "NS", decreasing = FALSE)], aes(b, -log10(p), colour = cl)) + geom_point(size = 0.3) +
  scale_colour_manual(values = c(`Tumor stronger` = cT, `Normal stronger` = cN, NS = "grey80"), name = NULL) +
  labs(x = "Tumor - normal, log odds", y = expression("-log10"~italic(P))) + theme_pub() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1)) + guides(colour = guide_legend(override.aes = list(size = 1.5)))

SP <- fread(file.path(D, "mapfilter_C_split.csv.gz"))
C3d <- rbindlist(lapply(c(1, 0.1, 0.01, 0.001), function(t) { z <- SP[p < t]; gs <- split(seq_len(nrow(z)), z$eg)
  bb <- replicate(2000, { i <- unlist(gs[sample(length(gs), replace = TRUE)], use.names = FALSE); mean(z$b[i] > 0) })
  data.table(lab = if (t == 1) "All" else as.character(t), n = nrow(z), f = mean(z$b > 0), lo = quantile(bb, .025), hi = quantile(bb, .975)) }))
C3d[, lab := factor(lab, levels = lab)]
note("Fig3C: fraction tumor-stronger by held-out P: %s", C3d[, paste(sprintf("%s %.3f [%.3f, %.3f] (n=%d)", lab, f, lo, hi, n), collapse = "; ")])
C3 <- ggplot(C3d, aes(lab, f)) + geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) + geom_pointrange(aes(ymin = lo, ymax = hi), colour = cT, size = 0.25, linewidth = 0.4) +
  coord_cartesian(ylim = c(0.45, 0.9)) + labs(x = expression("Held-out"~italic(P)~"threshold"), y = "Fraction tumor stronger") + theme_pub()

Y <- fread(file.path(V, "delta_ase", "dna_check_pairs_mapf.csv.gz"))[, grp := factor(grp, levels = c("tumor_stronger", "normal_stronger", "stable_control"))]
YL <- melt(Y[, .(grp, DNA = dnaShift, RNA = rnaShift)], id.vars = "grp")
note("Fig3D: SNP-patient pairs per group: %s", Y[, .N, by = grp][, paste(grp, N, collapse = "; ")])
D3 <- ggplot(YL, aes(grp, value, fill = variable)) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.2, width = 0.6, linewidth = 0.3, position = position_dodge(0.7)) + scale_fill_manual(values = c(DNA = "grey75", RNA = "#F4A582"), name = NULL) +
  coord_cartesian(ylim = c(-0.5, 0.5)) + scale_x_discrete(labels = c("Tumor\nstronger", "Normal\nstronger", "Stable")) +
  labs(x = NULL, y = "Tumor - normal shift") + theme_pub() + theme(legend.position = "top", legend.margin = margin(0, 0, -2, 0, "mm"))

M <- rbindlist(lapply(list.files(file.path(V, "delta_mbased_phased", "one"), full.names = TRUE), fread))[is.finite(maf_T) & is.finite(maf_N)]
MP <- M[, .(T = mean(maf_T), N = mean(maf_N)), by = pid]
note("Fig3E: MBASED phased one-sample, %d patients, %s patient-genes; MAF tumor %.4f vs normal %.4f; %.0f%% patients higher in tumor; Wilcoxon P=%.1g",
     nrow(MP), format(nrow(M), big.mark = ","), M[, mean(maf_T)], M[, mean(maf_N)], 100 * MP[, mean(T > N)], wilcox.test(MP$T, MP$N, paired = TRUE)$p.value)
limE <- range(c(MP$N, MP$T)); limE <- limE + c(-1, 1) * 0.05 * diff(limE)
E3 <- ggplot(MP, aes(N, T)) + geom_abline(linetype = 2, colour = "grey50", linewidth = 0.3) + geom_point(size = 0.8, colour = PAL[["dark"]], alpha = 0.8) +
  coord_equal(xlim = limE, ylim = limE, expand = FALSE) + labs(x = "Normal major-allele frequency", y = "Tumor major-allele frequency") + theme_pub()

T2 <- merge(rbindlist(lapply(list.files(file.path(V, "delta_mbased_phased", "two"), full.names = TRUE), fread)), M, by = c("pid", "eg"))
F3d <- rbindlist(lapply(c(0.05, 0.01, 0.001), function(t) { s <- T2[p_Tfirst < t]; pp <- s[, .(t = sum(maf_T > maf_N), n = sum(maf_T < maf_N)), by = pid]
  data.table(lab = as.character(t), n = nrow(s), f = s[, mean(maf_T > maf_N)], pT = pp[, sum(t > n)], pN = pp[, sum(n > t)], P = binom.test(pp[, sum(t > n)], pp[, sum(t != n)])$p.value) }))
F3d[, lab := factor(lab, levels = lab)]
note("Fig3F: MBASED phased two-sample, differential ASE: %s", F3d[, paste(sprintf("P<%s n=%d, tumor more imbalanced %.1f%%, patients %d vs %d (sign test P=%.1g)", lab, n, 100 * f, pT, pN, P), collapse = "; ")])
F3 <- ggplot(F3d, aes(lab, f)) + geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) + geom_col(fill = cT, width = 0.6) +
  coord_cartesian(ylim = c(0.4, 0.7)) + labs(x = expression("Differential ASE"~italic(P)~"threshold"), y = "Fraction tumor more imbalanced") + theme_pub()

g3 <- (A3 | B3 | C3) / (D3 | E3 | F3) + plot_annotation(tag_levels = "A") & tags
save_panel(g3, file.path(O, "Fig3.pdf"), 178, 118)

lg <- readLines(file.path(V, "34c_delta_eqtl_dedup_lmm_clean.log"))
lmm <- function(lab, label, grp) { l <- grep(sprintf("^\\[%s_dedup\\] mixed model", lab), lg, value = TRUE)
  x <- as.numeric(regmatches(l, gregexpr("[-+]?[0-9]*\\.?[0-9]+(e-?[0-9]+)?", l))[[1]])
  m <- x[1]; d <- x[2]; lo <- x[3]; hi <- x[4]
  data.table(label = label, group = grp, ratio = (m + d) / m, lo = (m + lo) / m, hi = (m + hi) / m, src = l) }
A4d <- rbind(lmm("gtex_prostate", "GTEx prostate lead", "Primary"), lmm("gtex_and_normal", "GTEx and CPGEA normal", "Primary"),
             lmm("normal_map", "CPGEA normal map", "Same-cohort"), lmm("both_maps", "CPGEA both maps", "Same-cohort"), lmm("tumour_map", "CPGEA tumor map", "Same-cohort"))
for (i in seq_len(nrow(A4d))) note("Fig4A %s: %.3f [%.3f, %.3f]  <- %s", A4d$label[i], A4d$ratio[i], A4d$lo[i], A4d$hi[i], A4d$src[i])
fp <- function(d, col) ggplot(d, aes(ratio, factor(label, levels = rev(unique(label))), colour = .data[[col]])) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.25, linewidth = 0.4) + labs(x = "Tumor / normal", y = NULL) + theme_pub()
A4 <- fp(A4d, "group") + scale_colour_manual(values = c(Primary = "black", `Same-cohort` = "grey65"), name = NULL) + theme(legend.position = "bottom", legend.margin = margin(-2, 0, 0, 0, "mm"))

H <- fread(file.path(D, "h3_trend_splits.csv"))
HL <- melt(H[, .(split = .I, Weak = r1, Medium = r2, Strong = r3)], id.vars = "split", variable.name = "base", value.name = "ratio")
HM <- HL[, .(ratio = median(ratio)), by = base]
note("Fig4B: %d splits; median ratio weak/medium/strong %s; weak-strong %.3f [%.3f, %.3f], P=%.3g; log-ratio slope %.3f [%.3f, %.3f]; absolute weak-strong P=%.3g",
     nrow(H), paste(sprintf("%.3f", HM$ratio), collapse = "/"), median(H$rdiff), median(H$rdiff_lo), median(H$rdiff_hi), median(H$p_rdiff),
     median(H$slope), median(H$slope_lo), median(H$slope_hi), median(H$p_adiff))
B4 <- ggplot(HL, aes(base, ratio)) + geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3) + geom_line(aes(group = split), colour = "grey85", linewidth = 0.3) +
  geom_point(colour = "grey60", size = 0.6) + geom_point(data = HM, colour = cT, size = 2.2, shape = 18) +
  scale_x_discrete(labels = c("Weak", "Mid", "Strong")) + labs(x = "Normal baseline", y = "Tumor / normal") + theme_pub()

DC <- fread(file.path(D, "delta_ase_decomposition.csv"))
C4d <- DC[stratum %in% c("GTEx: GTEx eGene", "GTEx: GTEx tested, not eGene"), .(label = fifelse(stratum == "GTEx: GTEx eGene", "GTEx eGene", "Not eGene"), group = "g", ratio, lo, hi)]
note("Fig4C: %s", C4d[, paste(sprintf("%s %.3f [%.3f, %.3f]", label, ratio, lo, hi), collapse = "; ")])
C4 <- fp(C4d, "group") + scale_colour_manual(values = c(g = PAL[["dark"]]), guide = "none")

H4 <- fread(file.path(D, "h4_splits.csv"))[base == "all", .(ratio = median(ratio), lo = median(lo), hi = median(hi), n = .N), by = orient]
H4[, label := fifelse(grepl("tumour", orient), "Tumor-preferred allele", "Normal-preferred allele")]
note("Fig4D: %s (median over %d splits)", H4[, paste(sprintf("%s %.3f [%.3f, %.3f]", label, ratio, lo, hi), collapse = "; ")], H4$n[1])
D4 <- fp(H4, "label") + scale_colour_manual(values = c(`Tumor-preferred allele` = cT, `Normal-preferred allele` = cN), guide = "none")

S40 <- fread(file.path(D, "subset_ratios_40.csv")); S41 <- fread(file.path(D, "subset_ratios_41.csv"))
E4d <- rbind(S40[baseline == "low" & subset == "all", .(label = sprintf("All patients (%d)", patients), ratio, lo, hi)],
             S40[baseline == "low" & subset == "normal epi >= 0.6", .(label = sprintf("Epithelium-rich normal (%d)", patients), ratio, lo, hi)],
             S40[baseline == "low" & subset == "smallest |delta epi| tertile", .(label = sprintf("Similar epithelial fraction (%d)", patients), ratio, lo, hi)],
             S41[baseline == "low" & subset == "smallest |delta basal share| tertile", .(label = sprintf("Similar basal share (%d)", patients), ratio, lo, hi)],
             S41[baseline == "low" & subset == "smallest |delta tumour-derived share| tertile", .(label = sprintf("Similar tumor-state share (%d)", patients), ratio, lo, hi)])[, group := "g"]
note("Fig4E (weak-baseline sites): %s", E4d[, paste(sprintf("%s %.3f [%.3f, %.3f]", label, ratio, lo, hi), collapse = "; ")])
E4 <- fp(E4d, "group") + scale_colour_manual(values = c(g = PAL[["dark"]]), guide = "none")

NL <- fread(file.path(V, "cistrans_slope_26e/tn_slope_strict_null.csv"))[model == "r0"]
l26 <- grep("^r0 +delta", readLines(file.path(V, "26e_tn_cistrans_slope_ndepth_1000.log")), value = TRUE)
obsH5 <- as.numeric(sub(".*obs ([-+0-9.]+).*", "\\1", l26)); pH5 <- (sum(abs(NL$delta) >= abs(obsH5)) + 1) / (nrow(NL) + 1)
note("Fig4F: H5 tumor-normal slope difference %+.4f; %d whole-patient label permutations, null %+.4f +- %.4f; two-sided P=%.3f  <- %s", obsH5, nrow(NL), mean(NL$delta), sd(NL$delta), pH5, trimws(l26))
F4 <- ggplot(NL, aes(delta)) + geom_histogram(bins = 40, fill = "grey75", colour = NA) + geom_vline(xintercept = obsH5, colour = cT, linewidth = 0.5) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) + scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) + labs(x = "Slope difference, T - N", y = "Permutations") + theme_pub()

g4 <- ((A4 | B4 | C4) + plot_layout(widths = c(1.35, 1, 0.9))) / ((D4 | E4 | F4) + plot_layout(widths = c(0.75, 1.15, 1))) + plot_annotation(tag_levels = "A") & tags
save_panel(g4, file.path(O, "Fig4.pdf"), 178, 118)
writeLines(LS, file.path(O, "fig3_fig4_legend_stats.txt")); cat("DONE\n")
