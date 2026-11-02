#!/usr/bin/env Rscript
# _clean (2026-10-07): per-patient Delta-ASE / Delta-eQTL from the cleaned reruns (56; 27_clean, 34c_clean); output to delta_eqtl_clean.
# Patient-level tumour-vs-own-normal amplification of germline cis effects vs clinical features.
# Metrics: (1) Delta-ASE (27, mappability-filtered, LOO orientation): sT/sN and sT - sN ;
#          (2) Delta-eQTL (34, normal-map leads): mean (tumour - normal) signed aFC per patient.
# Endpoints: Gleason score (Spearman; partial, adjusted for purity), PSA (Spearman), pathological T3+ vs T2 (Wilcoxon),
#            BCR (Cox per SD: univariable and adjusted for Gleason + purity). Small cohort (~117, ~30 BCR events): exploratory.
suppressPackageStartupMessages({library(data.table); library(survival)})
V <- "/path/to/data_root/ASE/code/v8_round2"; msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
cl <- fread("/path/to/data_root/CPGEA/CPGEA_web/CPGEA_Clinical_Info.tsv")
cl <- cl[, .(pid = as.character(`Sequencing ID`), gs = suppressWarnings(as.numeric(GS)), psa = suppressWarnings(as.numeric(PSA)),
             pT3 = grepl("T3|T4", `Pathology T Stage`), bcr = as.integer(toupper(BCR) == "YES"),
             months = suppressWarnings(as.numeric(gsub("[^0-9.]", "", `Month to surgery/BCR`))))]
A <- fread(file.path(V, "delta_ase_clean/mapfilter_A_per_patient.csv"))[, .(pid = as.character(pid), ase_ratio = sT / sN, ase_diff = sT - sN, purity)]
E <- fread(file.path(V, "delta_eqtl_clean/delta_eqtl_patients_normal_map.csv"))[, .(pid = as.character(pid), eqtl_diff = d, k)]
D <- merge(merge(A, E, by = "pid", all = TRUE), cl, by = "pid", all.x = TRUE)
msg("patients: Delta-ASE %d, Delta-eQTL %d, both %d ; with Gleason %d ; BCR events %d", D[is.finite(ase_diff), .N], D[is.finite(eqtl_diff), .N],
    D[is.finite(ase_diff) & is.finite(eqtl_diff), .N], D[is.finite(gs), .N], D[is.finite(ase_diff), sum(bcr, na.rm = TRUE)])
msg("Delta-ASE vs Delta-eQTL per patient: Spearman %+.3f (P=%.2g)", cor(D$ase_diff, D$eqtl_diff, method = "s", use = "c"),
    cor.test(D$ase_diff, D$eqtl_diff, method = "s", exact = FALSE)$p.value)
pc <- function(x, y, z) { ok <- is.finite(x) & is.finite(y) & is.finite(z); rx <- resid(lm(rank(x[ok]) ~ z[ok])); ry <- resid(lm(rank(y[ok]) ~ z[ok]))
  ct <- cor.test(rx, ry); c(r = unname(ct$estimate), p = ct$p.value, n = sum(ok)) }
for (m in c("ase_ratio", "ase_diff", "eqtl_diff")) {
  x <- D[[m]]; msg("\n== %s", m)
  g <- cor.test(x, D$gs, method = "s", exact = FALSE); gp <- pc(x, D$gs, D$purity)
  msg("  Gleason: Spearman %+.3f (P=%.3g) ; partial | purity %+.3f (P=%.3g, n=%d)", g$estimate, g$p.value, gp["r"], gp["p"], gp["n"])
  ps <- cor.test(x, D$psa, method = "s", exact = FALSE); msg("  PSA: Spearman %+.3f (P=%.3g)", ps$estimate, ps$p.value)
  ok <- is.finite(x); msg("  pT3+ vs pT2: median %.4f vs %.4f, Wilcoxon P=%.3g", median(x[ok & D$pT3]), median(x[ok & !D$pT3]), wilcox.test(x[ok & D$pT3], x[ok & !D$pT3])$p.value)
  S <- D[is.finite(get(m)) & is.finite(months) & is.finite(bcr)][, z := as.numeric(scale(get(m)))]
  c1 <- summary(coxph(Surv(months, bcr) ~ z, data = S))$coef; c2 <- summary(coxph(Surv(months, bcr) ~ z + gs + purity, data = S))$coef
  msg("  BCR (n=%d, events %d): HR per SD %.2f (P=%.3g) ; adjusted for Gleason + purity %.2f (P=%.3g)", nrow(S), sum(S$bcr), exp(c1["z", 1]), c1["z", 5], exp(c2["z", 1]), c2["z", 5]) }
fwrite(D, file.path(V, "delta_eqtl_clean/delta_clinical_patient_table.csv")); msg("DONE")
