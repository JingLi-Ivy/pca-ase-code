#!/usr/bin/env Rscript
# Fig 5B panels: endogenous, phased, tumour-DNA-corrected KIF18B allelic output at rs9914422 heterozygotes (data: 08_kif18b_phasing / 08b).
# 2026-10-10 rewrite: all numbers computed here (no hard-coded values). Points and band are on the same DNA-corrected scale:
#   per-tumour corrected C fraction = plogis(logit(RNA C fraction) - logit(tumour DNA C fraction)), DNA fraction clamped to [0.05, 0.95] as in 08b;
#   band = pooled beta-binomial intercept with DNA offset (same likelihood as 08b). aod::betabin ignores offset() in its formula, so the
#   hand-written likelihood is cross-checked against aod::betabin on the no-offset fit (must agree).
# A corrected C fraction vs FOXM1 (per-tumour Clopper-Pearson interval shifted by the same offset; size = phased reads)
# B RNA vs tumour-DNA C fraction (uncorrected)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(aod)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
D <- file.path(V, "kif18b_phase"); setwd(D)
P <- fread("phased_patient_level.csv")[rnaC + rnaG > 0 & is.finite(FOXM1)]
P[, n := rnaC + rnaG][, off := qlogis(pmin(pmax(dna_frac_C, 0.05), 0.95))][, FOXM1_z := as.numeric(scale(log2(FOXM1)))]
P[, state := fifelse(tcn == 2 & lcn == 1, "Balanced diploid", "SCNA / LOH")]

bb <- function(d, x = NULL) {
  y <- d$rnaC; n <- d$n; X <- cbind(1, if (!is.null(x)) d[[x]])
  nll <- function(th) { k <- ncol(X); mu <- plogis(X %*% th[1:k] + d$off); rho <- plogis(th[k + 1])
    a <- mu * (1 - rho) / rho; b <- (1 - mu) * (1 - rho) / rho; -sum(lchoose(n, y) + lbeta(y + a, n - y + b) - lbeta(a, b)) }
  f <- optim(c(rep(0, ncol(X)), -2), nll, hessian = TRUE, method = "BFGS"); se <- sqrt(diag(solve(f$hessian)))
  list(est = f$par, se = se, P = 2 * pnorm(-abs(f$par / se))) }
m0 <- bb(P); m1 <- bb(P, "FOXM1_z")
est <- plogis(m0$est[1]); ci <- plogis(m0$est[1] + c(-1, 1) * qnorm(0.975) * m0$se[1])
bal <- P[state == "Balanced diploid"]; m1b <- bb(bal, "FOXM1_z")
mz <- bb(copy(P)[, off := 0]); ab <- betabin(cbind(rnaC, rnaG) ~ 1, ~ 1, data = as.data.frame(P))
cat(sprintf("check (no offset): optim intercept %.4f (SE %.4f) vs aod::betabin %.4f (SE %.4f)\n", mz$est[1], mz$se[1], coef(ab)[1], sqrt(vcov(ab)[1, 1])))
stopifnot(abs(mz$est[1] - coef(ab)[1]) < 1e-3)

P[, c("lo", "hi") := as.list(binom.test(rnaC, n)$conf.int), by = pid]
adj <- function(p) plogis(qlogis(pmin(pmax(p, 0.005), 0.995)) - P$off)
P[, `:=`(fc = adj((rnaC + 0.5) / (n + 1)), flo = adj(lo), fhi = adj(hi))]
leg <- c(sprintf("tumours %d (balanced diploid %d), phased RNA reads %d (C %d, G %d; raw C fraction %.3f)", nrow(P), nrow(bal), sum(P$n), sum(P$rnaC), sum(P$rnaG), sum(P$rnaC) / sum(P$n)),
         sprintf("DNA-corrected C fraction (beta-binomial intercept) %.3f [%.3f, %.3f], P = %.2f", est, ci[1], ci[2], m0$P[1]),
         sprintf("FOXM1 slope (log2 FOXM1, z) %+.3f (SE %.3f), P = %.2f ; balanced diploid only (n = %d) P = %.2f", m1$est[2], m1$se[2], m1$P[2], nrow(bal), m1b$P[2]))
writeLines(leg, "Fig5_KIF18B_phased_legend_stats.txt"); cat(leg, sep = "\n")

cols <- c(`Balanced diploid` = PAL[["dark"]], `SCNA / LOH` = PAL[["accent"]])
th <- theme_pub()
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
A <- ggplot(P, aes(FOXM1, fc)) +
  annotate("rect", xmin = min(P$FOXM1) / 1.15, xmax = max(P$FOXM1) * 1.15, ymin = ci[1], ymax = ci[2], fill = "grey88") +
  geom_hline(yintercept = est, colour = PAL[["dark"]], linewidth = 0.4) + geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) +
  geom_errorbar(aes(ymin = flo, ymax = fhi, colour = state), width = 0, alpha = 0.4, linewidth = 0.3) + geom_point(aes(size = n, colour = state)) +
  scale_x_log10() + scale_y_continuous(limits = c(0, 1)) + scale_size_area(max_size = 2.5, breaks = c(10, 25, 50), name = "Phased reads") + scale_colour_manual(values = cols, name = NULL) +
  labs(x = "FOXM1 expression (FPKM)", y = "C-haplotype fraction\n(DNA-corrected)") + th + theme(legend.position = "right", legend.box.spacing = unit(1, "mm"))
B <- ggplot(P, aes(dna_frac_C, rnaC / n)) + geom_abline(linetype = 2, linewidth = 0.3) + geom_point(aes(size = n, colour = state)) +
  scale_size_area(max_size = 2.5, guide = "none") + scale_colour_manual(values = cols, guide = "none") + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Tumor DNA C fraction", y = "RNA C fraction") + th
g <- (A + B + plot_layout(widths = c(2, 1))) + plot_annotation(tag_levels = "A") & tags
save_panel(g, file.path(D, "Fig5_KIF18B_phased.pdf"), 178, 62)
save_panel(A, file.path(D, "Fig5B_KIF18B_phased_single.pdf"), 85, 58)
cat("DONE\n")
