#!/usr/bin/env Rscript
# 45c: pre-specified clinical test of the per-patient "awakened" gain (weak-baseline Delta-ASE gain, from 40's
# dilution_lowbaseline_patients.csv), saved as a script (the 2026-10-06 numbers were computed ad hoc; checklist line 528).
# Primary: Spearman(gain, Gleason) and partial on purity. Exploratory: BCR Cox per SD of rank-inverse-normal gain (raw and adjusted
# for Gleason + purity), median-split log-rank. Run on the original (delta_ase) and the cleaned (delta_ase_clean) per-patient tables.
suppressPackageStartupMessages({library(data.table); library(survival)})
V <- "/path/to/data_root/ASE/code/v8_round2"
cl <- fread("/path/to/data_root/CPGEA/CPGEA_web/CPGEA_Clinical_Info.tsv")[, .(pid = as.character(`Sequencing ID`), gs = suppressWarnings(as.numeric(GS)),
  bcr = fifelse(BCR == "YES", 1L, fifelse(BCR == "NO", 0L, NA_integer_)), months = suppressWarnings(as.numeric(sub("m$", "", `Month to surgery/BCR`))))]
pur <- fread(file.path(V, "delta_eqtl_clean/delta_clinical_patient_table.csv"))[, .(pid = as.character(pid), purity)]
int <- function(x) qnorm((rank(x) - 0.5) / length(x))
for (src in c("delta_ase", "delta_ase_clean")) {
  d <- merge(merge(fread(file.path(V, src, "dilution_lowbaseline_patients.csv"))[, .(pid = as.character(pid), gain)], cl, by = "pid"), pur, by = "pid")
  g <- d[is.finite(gs)]; s <- cor.test(g$gain, g$gs, method = "spearman", exact = FALSE)
  r <- residuals(lm(rank(gain) ~ purity, g)); q <- residuals(lm(rank(gs) ~ purity, g)); pc <- cor.test(r, q)
  b <- d[is.finite(bcr) & is.finite(months) & months > 0][, z := int(gain)]
  m0 <- summary(coxph(Surv(months, bcr) ~ z, b)); m1 <- summary(coxph(Surv(months, bcr) ~ z + gs + purity, b))
  b[, hi := gain >= median(gain)]; lr <- survdiff(Surv(months, bcr) ~ hi, b)
  cat(sprintf("[%s] n=%d | Gleason Spearman rho %+.3f (P=%.2f) ; partial|purity %+.3f (P=%.2f) | BCR n=%d events %d: HR/SD (rank-INT) %.2f [%.2f, %.2f] P=%.2f ; adjusted %.2f P=%.2f ; median split events %d vs %d, log-rank P=%.2f\n",
    src, nrow(g), s$estimate, s$p.value, pc$estimate, pc$p.value, nrow(b), sum(b$bcr), m0$conf.int[1, 1], m0$conf.int[1, 3], m0$conf.int[1, 4], m0$coef[1, 5],
    m1$conf.int[1, 1], m1$coef[1, 5], b[hi == TRUE, sum(bcr)], b[hi == FALSE, sum(bcr)], 1 - pchisq(lr$chisq, 1)))
}
