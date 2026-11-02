#!/usr/bin/env Rscript
# Step 1: are the Delta-ASE "tumor-stronger" SNPs explained by (subclonal) copy-number imbalance in tumor DNA?
# DNA allele counts from FACETS snp-pileup (File1 = normal, File2 = tumor), read-only inputs.
# For each patient x SNP pair, orient DNA to the RNA-defined high allele (dir from per-SNP test) and compare
# tumor-DNA vs normal-DNA oriented deviation; control = equal number of random stable SNPs.
suppressPackageStartupMessages(library(data.table)); setDTthreads(8); set.seed(5)
O <- "/path/to/data_root/ASE/code/v8_round2/delta_ase"; msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
IN <- Sys.getenv("DASE_IN", "delta_ase_per_snp.csv.gz"); TAG <- Sys.getenv("DASE_TAG", "")
R <- fread(file.path(O, IN))
B <- fread(file.path(O, "paired_balanced_sites.csv.gz")); B[, snp := paste(chrom, pos, ref, alt, sep = ":")]; B <- unique(B, by = c("pid","snp"))
sig <- if ("sig" %in% names(R)) R[sig == TRUE] else R[fdr < 0.05]; nsig <- nrow(sig); ctl <- R[fdr > 0.5][sample(.N, min(.N, 3 * nsig))]
sel <- rbind(sig[, .(snp, grp = ifelse(b > 0, "tumor_stronger", "normal_stronger"), dirR = sign(avg), b)],
             ctl[, .(snp, grp = "stable_control", dirR = sign(avg), b)])
X <- merge(B[, .(pid, snp, chrom, pos, ref, alt, rnT, anT, rnN, anN, purity)], sel, by = "snp")
msg("pairs: %s (%s)", format(nrow(X), big.mark = ","), paste(capture.output(print(X[, .N, by = grp])), collapse = " | "))
PD <- paste0("dna_pileup", TAG); dir.create(file.path(O, PD), showWarnings = FALSE)
for (p in unique(X$pid)) {
  pos <- unique(X[pid == p, .(chr = paste0("chr", chrom), pos)]); pf <- file.path(O, PD, paste0("pos_", p, ".txt"))
  fwrite(pos, pf, sep = ",", col.names = FALSE)
  out <- file.path(O, PD, paste0("T", p, ".csv"))
  if (!file.exists(out)) system(sprintf("zcat /path/to/data_root/variants/batch_43/facets/pileup/T%s.pileup.csv.gz | awk -F, 'NR==FNR{k[$1\",\"$2]=1;next} ($1\",\"$2) in k' %s - > %s", p, pf, out)) }
D <- rbindlist(lapply(unique(X$pid), function(p) { f <- file.path(O, PD, paste0("T", p, ".csv"))
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  d <- fread(f, header = FALSE, col.names = c("chr","pos","ref","alt","NR","NA_","NE","ND","TR","TA","TE","TD")); d[, pid := as.character(p)]; d }))
D[, chrom := sub("^chr", "", chr)]
X[, pid := as.character(pid)]; X[, chrom := as.character(chrom)]
Y <- merge(X, D[, .(pid, chrom, pos, dref = ref, dalt = alt, NR, NA_, TR, TA)], by = c("pid","chrom","pos"), all.x = TRUE)
msg("pairs with DNA pileup: %s of %s (%.1f%%)", format(Y[is.finite(TR), .N], big.mark = ","), format(nrow(Y), big.mark = ","), 100 * Y[, mean(is.finite(TR))])
Y <- Y[is.finite(TR) & (TR + TA) >= 10 & (NR + NA_) >= 10]
Y <- Y[(dref == ref & dalt == alt) | dref == "."]
Y[, `:=`(fTd = TR / (TR + TA), fNd = NR / (NR + NA_))]
Y[, `:=`(dnaT = dirR * (fTd - 0.5), dnaN = dirR * (fNd - 0.5), rnaT = dirR * (rnT / (rnT + anT) - 0.5), rnaN = dirR * (rnN / (rnN + anN) - 0.5))]
Y[, `:=`(dnaShift = dnaT - dnaN, rnaShift = rnaT - rnaN)]
Y[, p_dna := mapply(function(k, n, f) binom.test(k, n, min(max(f, 0.05), 0.95))$p.value, TR, TR + TA, fNd)]
msg("\nDNA vs RNA shift (oriented to RNA high allele), median [IQR]:")
print(Y[, .(pairs = .N, snps = uniqueN(snp), dna_shift = median(dnaShift), dna_q1 = quantile(dnaShift, .25), dna_q3 = quantile(dnaShift, .75),
            rna_shift = median(rnaShift), frac_dna_imbalanced_p01 = mean(p_dna < 0.01)), by = grp])
for (g in c("tumor_stronger","normal_stronger")) {
  a <- Y[grp == g, dnaShift]; c0 <- Y[grp == "stable_control", dnaShift]
  msg("%s vs stable control: tumor-DNA shift mean %+.4f vs %+.4f, Wilcoxon P=%.3g", g, mean(a), mean(c0), wilcox.test(a, c0)$p.value)
  t <- cor.test(Y[grp == g, dnaShift], Y[grp == g, rnaShift], method = "s", exact = FALSE)
  msg("%s: Spearman(DNA shift, RNA shift) %+.3f P=%.3g", g, t$estimate, t$p.value) }
agg <- Y[grp == "tumor_stronger", .(pairs = .N, dna_imb = sum(p_dna < 0.01), rna_all = mean(rnaShift), rna_dna_balanced = mean(rnaShift[p_dna >= 0.05]),
                                    n_bal = sum(p_dna >= 0.05), dna_mean = mean(dnaShift)), by = snp]
msg("\ntumor-stronger SNPs with DNA data: %d ; median share of pairs with DNA imbalance (P<0.01): %.2f", nrow(agg), median(agg$dna_imb / agg$pairs))
msg("RNA tumor-normal shift: all pairs %+.4f ; DNA-balanced pairs only %+.4f (SNPs with >=3 balanced pairs: %d)",
    agg[, mean(rna_all)], agg[n_bal >= 3, mean(rna_dna_balanced)], agg[n_bal >= 3, .N])
fwrite(Y, file.path(O, paste0("dna_check_pairs", TAG, ".csv.gz"))); fwrite(agg, file.path(O, paste0("dna_check_tumor_stronger_snps", TAG, ".csv")))
msg("DONE")
