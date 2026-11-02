#!/usr/bin/env Rscript
# _clean (2026-10-06): input delta_ase_clean/paired_balanced_sites.csv.gz (56: overflow rows and wrap-risk genes removed), weights capped at 100 reads; outputs to delta_ase_clean / delta_eqtl_clean.
# Validation of "tumour-awakened germline cis effects" (38: Delta-ASE gain concentrated at sites with little imbalance in normal tissue
# and at genes without a GTEx eQTL; direction shared across patients). Same filtered sites as 27/38.
# (1) CPGEA eQTL status of the host gene: tumour-only eGene / normal-only / both / neither -> patient-level T/N ratio (LOO, as 38).
#     Note: CPGEA maps come from the same patients' total expression (consistency check, not independent evidence).
# (2a) split-sample: per SNP (>= 8 patients) random halves; orientation and baseline (|normal deviation|) from half 1 only,
#      tumour and normal deviations measured in half 2 only.
# (2b) cross-tissue orientation in the split design: orient by half-1 tumour only, or by half-1 normal only; tumour-awakened effects
#      predict T2 >> N2 when oriented by tumour, and T2 ~ N2 when oriented by normal.
# (2c) technical check: at low-baseline sites, is the tumour-favoured allele (half-1 tumour orientation) the reference allele more
#      often than 50%? (reference/mapping bias specific to tumour libraries would push it above 50%.)
suppressPackageStartupMessages(library(data.table)); setDTthreads(8); set.seed(39)
V <- "/path/to/data_root/ASE/code/v8_round2"; D <- file.path(V, "delta_ase_clean"); F <- file.path(V, "mapfilter")
QD <- "/path/to/data_root/CPGEA/QTLs/CPGEA_eQTL_and_sQTL/"; msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
B <- fread(file.path(D, "paired_balanced_sites.csv.gz"))
B[, `:=`(snp = paste(chrom, pos, ref, alt, sep = ":"), chrom = as.character(chrom), eg = sub("\\..*", "", eg))]
sd <- fread(file.path(F, "genomicSuperDups.txt.gz"), select = 2:4, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
bl <- fread(file.path(F, "hg38-blacklist.v2.bed.gz"), select = 1:3, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
ov <- function(reg) { q <- unique(B[, .(chrom, s = pos, e = pos)]); r <- unique(reg[, .(chrom, s = s + 1L, e)]); setkey(r, chrom, s, e)
  foverlaps(q, r, nomatch = 0L)[, unique(paste(chrom, i.s))] }
bad <- c(ov(sd), ov(bl))
gtf <- fread(cmd = "zcat /path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz | awk -F'\\t' '$3==\"gene\"' | cut -f9", header = FALSE, sep = "\t")
pseudo <- gtf[grepl('gene_biotype "[a-z_]*pseudogene"', V1), sub('.*gene_id "([^"]+)".*', "\\1", V1)]
mono <- unique(B, by = c("pid", "snp"))[, .(devN = abs(sum(rnN) / sum(rnN + anN) - 0.5)), by = snp][devN > 0.35, snp]
B <- B[!(paste(chrom, pos) %chin% bad | (chrom == "6" & pos >= 28510120 & pos <= 33480577) | eg %chin% pseudo | snp %chin% mono)]
B <- unique(B, by = c("pid", "snp")); B[, `:=`(dT = rnT + anT, dN = rnN + anN)]

B[, npat := .N, by = snp]
one_split <- function(seed) { set.seed(seed)
  S <- B[npat >= 8][, half := sample(rep(1:2, length.out = .N)), by = snp]
  H1 <- S[half == 1, .(dT1 = sum(rnT) / sum(dT) - 0.5, dN1 = sum(rnN) / sum(dN) - 0.5), by = snp]
  H2 <- merge(S[half == 2], H1, by = "snp")[, o := 0.5 * (dT1 + dN1)][abs(o) >= 0.05]
  H2[, base := cut(abs(dN1), quantile(abs(dN1), 0:3 / 3), include.lowest = TRUE, labels = 1:3)]
  H2[, `:=`(sT = sign(o) * (rnT / dT - 0.5), sN = sign(o) * (rnN / dN - 0.5), w = pmin(dN, dT, 100))]
  P <- H2[, .(sT = sum(w * sT), sN = sum(w * sN), W = sum(w)), by = .(pid, base)]
  stat <- function(P) { a <- P[, .(T = sum(sT) / sum(W), N = sum(sN) / sum(W)), by = base][order(base)]
    r <- a$T / a$N; c(rdiff = r[1] - r[3], slope = unname(coef(lm(log(r) ~ c(1, 2, 3)))[2]), adiff = (a$T[1] - a$N[1]) - (a$T[3] - a$N[3]), r1 = r[1], r2 = r[2], r3 = r[3],
      g1 = a$T[1] - a$N[1], g2 = a$T[2] - a$N[2], g3 = a$T[3] - a$N[3]) }
  ob <- stat(P); pids <- unique(P$pid)
  bt <- t(replicate(1000, { s <- sample(pids, replace = TRUE); stat(rbindlist(lapply(seq_along(s), function(i) P[pid == s[i]]))) }))
  c(ob, rdiff_lo = quantile(bt[, "rdiff"], .025, names = FALSE), rdiff_hi = quantile(bt[, "rdiff"], .975, names = FALSE), p_rdiff = mean(bt[, "rdiff"] <= 0),
    slope_lo = quantile(bt[, "slope"], .025, names = FALSE), slope_hi = quantile(bt[, "slope"], .975, names = FALSE), p_slope = mean(bt[, "slope"] >= 0),
    adiff_lo = quantile(bt[, "adiff"], .025, names = FALSE), adiff_hi = quantile(bt[, "adiff"], .975, names = FALSE), p_adiff = 2 * min(mean(bt[, "adiff"] <= 0), mean(bt[, "adiff"] >= 0))) }
R <- do.call(rbind, parallel::mclapply(1:20, function(k) one_split(3900 + k), mc.cores = 10))
fwrite(as.data.table(R), file.path(D, "h3_trend_splits.csv"))
m <- colMeans(R); md <- apply(R, 2, median)
msg("H3 trend test, 20 random splits (medians across splits of per-split estimates / bootstrap CIs / one-sided P):")
msg("  ratio by baseline tertile weak / mid / strong: %.3f / %.3f / %.3f", md["r1"], md["r2"], md["r3"])
msg("  PRIMARY  ratio(weak) - ratio(strong) = %+.3f [%+.3f, %+.3f], P = %.3g (range of P across splits %.3g-%.3g)", md["rdiff"], md["rdiff_lo"], md["rdiff_hi"], md["p_rdiff"], min(R[, "p_rdiff"]), max(R[, "p_rdiff"]))
msg("  trend    slope of log(ratio) per tertile = %+.3f [%+.3f, %+.3f], P = %.3g", md["slope"], md["slope_lo"], md["slope_hi"], md["p_slope"])
msg("  absolute gain (tumour - normal) weak / mid / strong: %+.4f / %+.4f / %+.4f ; weak - strong = %+.4f [%+.4f, %+.4f], two-sided P = %.3g",
    md["g1"], md["g2"], md["g3"], md["adiff"], md["adiff_lo"], md["adiff_hi"], md["p_adiff"])
msg("DONE")
