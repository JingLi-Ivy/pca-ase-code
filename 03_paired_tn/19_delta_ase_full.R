#!/usr/bin/env Rscript
# Delta-ASE, full version. Input: paired tumor-normal sites (tumor balanced diploid, >= 10 reads in both) from 18_delta_ase_pilot.R.
# Orientation and site selection never favour one tissue:
#  (A) patient-level aggregate: each patient's direction at a SNP comes from OTHER patients (leave-one-out),
#      using the equal-weight average of tumor and normal pooled deviations; selection |LOO dev| >= DEV_MIN in >= K others.
#  (B) per-SNP test: selection/orientation by the equal-weight T/N average (symmetric); paired quasi-binomial GLM
#      cbind(high, low) ~ patient + tissue ; tissue coefficient > 0 = stronger in tumor ("tumor-gained"), < 0 = "tumor-lost".
suppressPackageStartupMessages(library(data.table)); setDTthreads(16); set.seed(3)
O <- "/path/to/data_root/ASE/code/v8_round2/delta_ase"; msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
K <- 5L; DEV_MIN <- 0.05
B <- fread(file.path(O, "paired_balanced_sites.csv.gz"))
B[, snp := paste(chrom, pos, ref, alt, sep = ":")]
B <- unique(B, by = c("pid", "snp"))
msg("paired balanced patient-SNP observations: %s ; SNPs %s ; patients %d", format(nrow(B), big.mark = ","), format(uniqueN(B$snp), big.mark = ","), uniqueN(B$pid))
B[, npat := .N, by = snp]; B <- B[npat >= K + 1]
msg("SNPs heterozygous & covered in >= %d patients: %s (%s obs)", K + 1, format(uniqueN(B$snp), big.mark = ","), format(nrow(B), big.mark = ","))

B[, dT := rnT + anT][, dN := rnN + anN]
tot <- B[, .(RT = sum(rnT), DT = sum(dT), RN = sum(rnN), DN = sum(dN)), by = snp]
B <- merge(B, tot, by = "snp")
B[, loo := 0.5 * ((RT - rnT) / (DT - dT) + (RN - rnN) / (DN - dN)) - 0.5]
A <- B[abs(loo) >= DEV_MIN]
A[, `:=`(sN = sign(loo) * (rnN / dN - 0.5), sT = sign(loo) * (rnT / dT - 0.5))]
pat <- A[, .(sN = weighted.mean(sN, pmin(dN, dT)), sT = weighted.mean(sT, pmin(dN, dT)), n = .N, purity = purity[1]), by = pid]
bt <- replicate(2000, { p <- pat[sample(.N, replace = TRUE)]; c(p[, weighted.mean(sN, n)], p[, weighted.mean(sT, n)]) })
mN <- pat[, weighted.mean(sN, n)]; mT <- pat[, weighted.mean(sT, n)]
msg("\n(A) LOO-oriented sites: %s obs, %s SNPs, %d patients", format(nrow(A), big.mark = ","), format(uniqueN(A$snp), big.mark = ","), uniqueN(A$pid))
msg("(A) mean oriented deviation: normal %+.4f [%+.4f,%+.4f] ; tumor %+.4f [%+.4f,%+.4f] ; tumor/normal %.3f [%.3f,%.3f]",
    mN, quantile(bt[1, ], .025), quantile(bt[1, ], .975), mT, quantile(bt[2, ], .025), quantile(bt[2, ], .975),
    mT / mN, quantile(bt[2, ] / bt[1, ], .025), quantile(bt[2, ] / bt[1, ], .975))
pat[, pq := cut(purity, quantile(purity, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low","mid","high"))]
print(pat[, .(patients = .N, normal = weighted.mean(sN, n), tumor = weighted.mean(sT, n)), by = pq][order(pq)][, ratio := tumor / normal][])
ct <- cor.test(pat$sT / pat$sN, pat$purity, method = "s", exact = FALSE)
msg("(A) per-patient tumor/normal ratio vs purity: Spearman %+.3f P=%.3g (loss in tumor cells predicts a negative slope)", ct$estimate, ct$p.value)

S <- B[, .(avg = 0.5 * (sum(rnT) / sum(dT) + sum(rnN) / sum(dN)) - 0.5, np = .N), by = snp][abs(avg) >= DEV_MIN]
X <- B[snp %in% S$snp]; X <- merge(X, S[, .(snp, dir = sign(avg))], by = "snp")
X[, `:=`(hT = ifelse(dir > 0, rnT, anT), hN = ifelse(dir > 0, rnN, anN))]
L <- rbind(X[, .(snp, pid, tissue = 1L, h = hT, d = dT, purity)], X[, .(snp, pid, tissue = 0L, h = hN, d = dN, purity)])
msg("\n(B) SNPs tested: %s (%s patient-SNP pairs)", format(uniqueN(L$snp), big.mark = ","), format(nrow(X), big.mark = ","))
fit1 <- function(z) { f <- tryCatch(suppressWarnings(glm(cbind(h, d - h) ~ factor(pid) + tissue, family = quasibinomial, data = z)), error = function(e) NULL)
  if (is.null(f)) return(list(b = NA_real_, se = NA_real_, p = NA_real_)); s <- summary(f)$coef
  if (!"tissue" %in% rownames(s)) return(list(b = NA_real_, se = NA_real_, p = NA_real_)); list(b = s["tissue", 1], se = s["tissue", 2], p = s["tissue", 4]) }
setkey(L, snp)
R <- L[, fit1(.SD), by = snp]
R <- merge(R, S, by = "snp"); R <- merge(R, unique(B[, .(snp, chrom, pos, eg)], by = "snp"), by = "snp")
R[, fdr := p.adjust(p, "BH")]
R <- merge(R, X[, .(fT = sum(hT) / sum(dT) - 0.5, fN = sum(hN) / sum(dN) - 0.5), by = snp], by = "snp")
R[, class := fcase(fdr < 0.05 & b > 0 & fN > 0.02, "tumor-magnified", fdr < 0.05 & b > 0 & fN <= 0.02, "tumor-gained",
                   fdr < 0.05 & b < 0 & fT > 0.02, "normal-magnified", fdr < 0.05 & b < 0 & fT <= 0.02, "tumor-lost", default = "stable")]
msg("(B) tissue effect FDR<0.05: %d of %d SNPs (%.1f%%)", R[fdr < 0.05, .N], nrow(R), 100 * R[fdr < 0.05, .N] / nrow(R))
print(R[, .N, by = class][order(-N)][, pct := round(100 * N / sum(N), 1)][])
msg("(B) among significant: tumor stronger %d, normal stronger %d", R[fdr < 0.05 & b > 0, .N], R[fdr < 0.05 & b < 0, .N])
Y <- X[snp %in% R[fdr < 0.05, snp]]; Y <- merge(Y, R[, .(snp, b)], by = "snp")
Y[, dlt := (hT / dT - 0.5) - (hN / dN - 0.5)]
for (sg in c(1, -1)) { z <- Y[sign(b) == sg]; if (nrow(z) > 50) { t <- cor.test(z$dlt * sg, z$purity, method = "s", exact = FALSE)
  msg("(B) %s SNPs: |T-N difference| vs purity Spearman %+.3f P=%.3g (n=%d pairs)", ifelse(sg > 0, "tumor-stronger", "normal-stronger"), t$estimate, t$p.value, nrow(z)) } }
fwrite(R[order(p)], file.path(O, "delta_ase_per_snp.csv.gz")); fwrite(pat, file.path(O, "delta_ase_per_patient_LOO.csv"))
msg("DONE")
