#!/usr/bin/env Rscript
# Phase-resolved KIF18B allelic output: for each rs9914422-heterozygous tumor, assign RNA/DNA reads of each reported SNP to the C vs G haplotype
# A_reg = log2(RNA_C/RNA_G) - log2(DNA_C/DNA_G); formal inference uses a beta-binomial (overdispersed) model with the tumor DNA fraction as offset
suppressPackageStartupMessages(library(data.table))
setwd("/path/to/data_root/ASE/code/v8_round2/kif18b_phase"); A <- "/path/to/data_root/ASE/code/v8_round2/kif18b_audit"
RS <- 44947872
system("bcftools query -f '%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n' phased.vcf.gz > phased.tsv")
ph <- fread("phased.tsv", header = FALSE, sep = "\t")
setnames(ph, 1:3, c("pos","ref","alt"))
long <- melt(ph, id.vars = c("pos","ref","alt"), value.name = "x")[, c("sample","gt") := tstrsplit(x, "=", fixed = TRUE)]
long[, pid := gsub("\\D", "", sub("_.*", "", sample))][, c("h1","h2") := tstrsplit(gt, "|", fixed = TRUE)]
long[, `:=`(h1 = as.integer(h1), h2 = as.integer(h2))]
rs <- long[pos == RS, .(pid, rs_h1 = h1, rs_h2 = h2)]
het <- rs[rs_h1 != rs_h2]; cat("rs9914422 phased hets:", nrow(het), "\n")
sn <- merge(long[h1 != h2, .(pid, pos, ref, alt, h1, h2)], het, by = "pid")
sn[, C_on_alt := fifelse(rs_h1 == 1L, h1 == 1L, h2 == 1L)]
r <- fread(file.path(A, "kif18b_rna_snps.csv"))[gene_id == "ENSG00000186185"][, pid := gsub("\\D", "", sample_id)]
r <- r[, .(pid, pos = start, rref = ref, ralt = alt, Rr = ref.matches, Ra = alt.matches)]
p <- fread(file.path(A, "kif18b_dna_pileup.csv"), header = FALSE,
           col.names = c("sample","chr","pos","ref","alt","NR","NA_","NE","ND","TR","TA","TE","TD"))[, pid := gsub("\\D", "", sample)]
m <- merge(sn, r, by = c("pid","pos"))
m <- m[(rref == ref & ralt == alt)]
m <- merge(m, p[, .(pid, pos, TR, TA, NR, NA_)], by = c("pid","pos"), all.x = TRUE)
m[, `:=`(rnaC = fifelse(C_on_alt, Ra, Rr), rnaG = fifelse(C_on_alt, Rr, Ra),
         dnaC = fifelse(C_on_alt, TA, TR), dnaG = fifelse(C_on_alt, TR, TA))]
cat(sprintf("phased RNA SNP obs: %d across %d rs9914422-het tumors (SNPs %d); with DNA pileup: %d\n",
            nrow(m), uniqueN(m$pid), uniqueN(m$pos), m[!is.na(TR), .N]))
fwrite(m, "phased_snp_level.csv")
au <- fread(file.path(A, "kif18b_patient_audit_all.csv"))[, pid := as.character(pid)]
pt <- m[, .(n_snp = .N, rnaC = sum(rnaC), rnaG = sum(rnaG), dnaC = sum(dnaC, na.rm = TRUE), dnaG = sum(dnaG, na.rm = TRUE)), by = pid]
pt <- merge(pt, au[, .(pid, FOXM1, FOXM1_tertile, sig_ASE, tcn, lcn, purity, epithelial, prolif_score, KIF18B, KIF18B_DT)], by = "pid", all.x = TRUE)
pt[, `:=`(rna_frac_C = rnaC / (rnaC + rnaG), dna_frac_C = (dnaC + 0.5) / (dnaC + dnaG + 1))]
pt[, A_reg := log2((rnaC + 0.5) / (rnaG + 0.5)) - log2((dnaC + 0.5) / (dnaG + 0.5))]
fwrite(pt, "phased_patient_level.csv")
cat("\n===== Patient level (rs9914422 heterozygous, phased RNA reads) =====\n")
print(pt[order(-(rnaC + rnaG)), .(pid, n_snp, rnaC, rnaG, rna_frac_C = round(rna_frac_C, 2), dnaC, dnaG, dna_frac_C = round(dna_frac_C, 2), A_reg = round(A_reg, 2), FOXM1 = round(FOXM1, 1), tert = FOXM1_tertile, sig_ASE, tcn, lcn)])
tot <- pt[, .(C = sum(rnaC), G = sum(rnaG))]
cat(sprintf("\nTotal RNA reads: C haplotype %d vs G haplotype %d (C fraction %.3f; binomial P=%.3g)\n", tot$C, tot$G, tot$C / (tot$C + tot$G), binom.test(tot$C, tot$C + tot$G)$p.value))
cat(sprintf("Patients with C>G: %d/%d (sign test P=%.3g)\n", pt[rnaC > rnaG, .N], pt[rnaC != rnaG, .N], binom.test(pt[rnaC > rnaG, .N], pt[rnaC != rnaG, .N])$p.value))
bb <- function(d, x = NULL, offset = TRUE) {
  y <- d$rnaC; n <- d$rnaC + d$rnaG; off <- if (offset) qlogis(pmin(pmax(d$dna_frac_C, 0.05), 0.95)) else 0
  X <- cbind(1, if (!is.null(x)) d[[x]])
  nll <- function(th) { b <- th[1:ncol(X)]; rho <- plogis(th[ncol(X) + 1]); mu <- plogis(X %*% b + off)
    a <- mu * (1 - rho) / rho; bb_ <- (1 - mu) * (1 - rho) / rho
    -sum(lchoose(n, y) + lbeta(y + a, n - y + bb_) - lbeta(a, bb_)) }
  f <- optim(c(rep(0, ncol(X)), -2), nll, hessian = TRUE, method = "BFGS")
  se <- sqrt(diag(solve(f$hessian))); k <- ncol(X)
  data.table(term = c("intercept", if (k > 1) x, "logit_rho")[1:(k + 1)], est = f$par, se = se, z = f$par / se, P = 2 * pnorm(-abs(f$par / se)), nll = f$value) }
d <- pt[rnaC + rnaG >= 1 & is.finite(FOXM1)][, FOXM1_z := as.numeric(scale(log2(FOXM1)))]
cat(sprintf("\n===== beta-binomial (n=%d tumors, %d reads) =====\n", nrow(d), d[, sum(rnaC + rnaG)]))
cat("intercept only, DNA offset: intercept>0 = higher C-haplotype expression after DNA correction\n"); print(bb(d))
cat("FOXM1 (log2, z) slope, DNA offset: slope>0 = higher FOXM1, larger C/G ratio (amplification)\n"); print(bb(d, "FOXM1_z"))
cat("as above, without DNA offset\n"); print(bb(d, "FOXM1_z", offset = FALSE))
d2 <- d[tcn == 2 & lcn == 1]
if (nrow(d2) >= 8) { cat(sprintf("Locally balanced diploid only (tcn2 lcn1, n=%d)\n", nrow(d2))); print(bb(d2, "FOXM1_z")) }
cat(sprintf("\nSpearman(A_reg, FOXM1) = %.3f, P = %.3g (n=%d)\n", cor(d$A_reg, d$FOXM1, method = "s"), cor.test(d$A_reg, d$FOXM1, method = "s", exact = FALSE)$p.value, nrow(d)))
