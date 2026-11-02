#!/usr/bin/env Rscript
# _clean (2026-10-06): input delta_ase_clean/paired_balanced_sites.csv.gz (56: overflow rows and wrap-risk genes removed), weights capped at 100 reads; outputs to delta_ase_clean / delta_eqtl_clean.
# 34c: as 34b plus (i) read de-duplication: within each patient x gene, reporter SNPs are kept greedily by total depth so that kept
#      SNPs are >= 300 bp apart (a read pair is then counted at most once); (ii) formal test by linear mixed model
#      signed log2 aFC ~ tissue + (1|gene) + (1|patient) + (1|gene:patient) (lme4), with and without de-duplication.
#      GTEx lead set is used when delta_eqtl/Prostate.v8.egenes.txt.gz exists.
# 34b: as 34, adding an independent lead set from GTEx v8 Prostate (egenes, qval < 0.05; lead variant and ALT-allele slope),
#      which carries no winner's curse in either CPGEA tissue. GTEx ALT/slope enters as minor/beta, so the direction logic of 34 is reused.
# 34: within-patient Delta-eQTL. Set-up (filters, lead eQTLs, phased genotypes) copied from 31b.
# For each lead cis-eQTL and each eQTL-heterozygous patient with phased reporter SNPs, signed log2 aFC (eQTL direction) is computed in
# BOTH the tumour and the matched normal sample (same patient, same haplotypes, same reporter SNPs; tumour balanced-diploid sites).
# Lead maps: normal (winner's curse inflates normal effects -> biased towards "weaker in tumour"), tumour (opposite bias), and genes
# that are eQTLs in both maps (normal lead). Paired comparison per gene (Wilcoxon), global means with gene- and patient-level
# bootstrap, patient-level difference vs purity / epithelial fraction. Patient-level table written for the clinical analysis (35).
# Set-up copied from:
# cis x trans, Flynn et al. 2022 (PLoS Genet) main design, applied to CPGEA tumour and normal prostate.
#   For every lead cis-eQTL (CPGEA map, FDR < 0.05), among eQTL-heterozygous patients: phased haplotypic RNA counts summed over
#   the gene's reporter SNPs (Beagle-phased germline, 31_phase_genome.sh) -> log2 aFC of the haplotype carrying the eQTL ALT
#   allele, signed so that > 0 = imbalance in the eQTL direction. Then partial Spearman correlation of signed aFC with each TF
#   (log2 TPM), adjusting for epithelial fraction (+ purity in tumour). z > 0: more TF -> stronger eQTL effect.
#   Only eQTL heterozygotes are used (no het vs hom comparison), so unequal measurement reliability cannot create a signal.
# Sites: paired balanced-diploid sites (tumour tcn2/lcn1, >= 10 reads in both tissues) after the mappability filters of 27.
# Positive control: mean signed aFC > 0 (ASE direction agrees with the eQTL). Nulls: strict patient permutation of TF rows (100);
# expression-matched random non-TF genes (20 draws).
suppressPackageStartupMessages({library(data.table); library(parallel)}); setDTthreads(8); set.seed(31)
V <- "/path/to/data_root/ASE/code/v8_round2"; D <- file.path(V, "delta_ase_clean"); F <- file.path(V, "mapfilter"); PH <- file.path(V, "phase_genome")
O <- file.path(V, "delta_eqtl_clean"); dir.create(O, showWarnings = FALSE); QD <- "/path/to/data_root/CPGEA/QTLs/CPGEA_eQTL_and_sQTL/"
AD <- "/path/to/data_root/ASE/ADASTRA6/ADASTRA.v.6.1.Mabel/TF"
NPERM <- as.integer(Sys.getenv("NPERM", "100")); NDRAW <- as.integer(Sys.getenv("NDRAW", "20")); NC <- as.integer(Sys.getenv("NC", "16"))
MINHET <- 15L; MINREADS <- 10L
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

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
B[, pid := as.character(pid)]

lead <- function(tis) {
  q <- fread(paste0(QD, sprintf("CPGEA_cis_%s_eQTLs_anno.txt", if (tis == "T") "tumor" else "normal")))[FDR < 0.05]
  q[, eg := sub("\\..*", "", sub("\\|.*", "", gene))]; q <- q[order(p.value)][!duplicated(eg)][eg %chin% unique(B$eg)]
  q[, .(eg, chrom = sub("^chr", "", Chr), lpos = Pos, beta, minor = Minor.Allele, major = Major.Allele)] }
QT <- lead("T"); QN <- lead("N")
GF <- "/path/to/data_root/ASE/code/v8_round2/delta_eqtl/Prostate.v8.egenes.txt.gz"
if (file.exists(GF)) {
gt <- fread(GF)[qval < 0.05]
gt[, c("vchr", "vpos", "vref", "valt") := tstrsplit(variant_id, "_", keep = 1:4)]
QG <- gt[nchar(vref) == 1 & nchar(valt) == 1, .(eg = sub("\\..*", "", gene_id), chrom = sub("^chr", "", vchr), lpos = as.integer(vpos), beta = slope, minor = valt, major = vref)][eg %chin% unique(B$eg)]
msg("GTEx prostate eGenes (qval<0.05, SNV lead) with CPGEA reporter SNPs: %d", nrow(QG))
} else { QG <- QT[0]; msg("GTEx file not found: GTEx lead set skipped") }

need <- unique(rbind(QT[, .(chrom, pos = lpos)], QN[, .(chrom, pos = lpos)], QG[, .(chrom, pos = lpos)], B[eg %chin% c(QT$eg, QN$eg, QG$eg), .(chrom, pos)]))
pf <- file.path(O, "flynn_positions_34c.tsv"); fwrite(need[order(as.integer(chrom), pos)], pf, sep = "\t", col.names = FALSE)
GTL <- rbindlist(lapply(1:22, function(c) { f <- file.path(PH, sprintf("phased_chr%d.vcf.gz", c)); if (!file.exists(f)) return(NULL)
  smp <- system(sprintf("bcftools query -l %s", f), intern = TRUE)
  x <- fread(cmd = sprintf("bcftools query -R %s -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' %s", pf, f), header = FALSE, sep = "\t", colClasses = "character")
  if (!nrow(x)) return(NULL); setnames(x, c("chrom", "pos", "ref", "alt", gsub("\\D", "", smp)))
  melt(x, id.vars = c("chrom", "pos", "ref", "alt"), variable.name = "pid", value.name = "gt") }))
GTL[, `:=`(pos = as.integer(pos), pid = as.character(pid))][, `:=`(h1 = substr(gt, 1, 1), h2 = substr(gt, 3, 3))]
GTL <- GTL[gt %in% c("0|1", "1|0")]
msg("phased heterozygous genotypes at needed positions: %s", format(nrow(GTL), big.mark = ","))

pur <- fread("/path/to/data_root/variants/batch_43/facets/Table_S1_PRAD.txt")[, .(pid = as.character(SAMPLE_ID), purity = as.numeric(purity))]
epi <- fread(file.path(V, "decon_v2d/CPGEA_tumor_cell_fraction.tsv")); epi <- epi[, .(pid = sub("^T", "", get(names(epi)[1])), epiT = Epithelial_cells)]
paired <- function(Q, lab, DEDUP = FALSE) {
  Lg <- merge(Q, unique(GTL[, .(chrom, pos, ref, alt)]), by.x = c("chrom", "lpos"), by.y = c("chrom", "pos"))
  Lg <- Lg[(minor == alt & major == ref) | (minor == ref & major == alt)][, sgn := fifelse(minor == alt, sign(beta), -sign(beta))]
  Lh <- merge(Lg[, .(eg, chrom, lpos, sgn)], GTL[, .(chrom, lpos = pos, pid, hl1 = h1)], by = c("chrom", "lpos"))
  R <- merge(B[, .(pid, eg, chrom, pos, rnT, anT, rnN, anN)], GTL[, .(chrom, pos, pid, r1 = h1)], by = c("chrom", "pos", "pid"))
  R <- merge(R, Lh, by = c("eg", "pid", "chrom"))[, same := r1 == hl1]
  if (DEDUP) { R[, dep := rnT + anT + rnN + anN]; setorder(R, eg, pid, -dep)
    R <- R[, { kp <- logical(.N); kept <- integer(0); for (j in seq_len(.N)) if (!length(kept) || all(abs(pos[j] - kept) >= 300L)) { kp[j] <- TRUE; kept <- c(kept, pos[j]) }; .SD[kp] }, by = .(eg, pid)] }
  A <- R[, .(laT = sum(fifelse(same, anT, rnT)), lrT = sum(fifelse(same, rnT, anT)), laN = sum(fifelse(same, anN, rnN)), lrN = sum(fifelse(same, rnN, anN)), sgn = sgn[1]), by = .(eg, pid)]
  A <- A[laT + lrT >= MINREADS & laN + lrN >= MINREADS]
  A[, `:=`(aT = sgn * log2((laT + 0.5) / (lrT + 0.5)), aN = sgn * log2((laN + 0.5) / (lrN + 0.5)))][, d := aT - aN]
  A <- A[, if (.N >= 10) .SD, by = eg]
  G <- A[, .(n = .N, mT = mean(aT), mN = mean(aN), p = tryCatch(wilcox.test(aT, aN, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_)), by = eg][, fdr := p.adjust(p, "BH")]
  gb <- replicate(2000, { g <- G[sample(.N, replace = TRUE)]; c(mean(g$mT) / mean(g$mN), mean(g$mT - g$mN)) })
  Pp <- A[, .(d = mean(d), mT = mean(aT), mN = mean(aN), k = .N), by = pid]
  pb <- replicate(2000, { p <- Pp[sample(.N, replace = TRUE)]; weighted.mean(p$mT, p$k) / weighted.mean(p$mN, p$k) })
  msg("\n[%s] eQTL genes with >= 10 paired heterozygous patients: %d ; patient-gene pairs %s ; patients %d", lab, nrow(G), format(nrow(A), big.mark = ","), uniqueN(A$pid))
  msg("[%s] mean signed log2 aFC: tumour %+.4f, normal %+.4f ; tumour/normal %.3f (gene bootstrap [%.3f, %.3f] ; patient bootstrap [%.3f, %.3f]) ; mean difference %+.4f [%+.4f, %+.4f]",
      lab, mean(G$mT), mean(G$mN), mean(G$mT) / mean(G$mN), quantile(gb[1, ], .025), quantile(gb[1, ], .975), quantile(pb, .025), quantile(pb, .975),
      mean(G$mT - G$mN), quantile(gb[2, ], .025), quantile(gb[2, ], .975))
  bt <- binom.test(sum(G$mT > G$mN), nrow(G))
  msg("[%s] genes with stronger eQTL effect in tumour %.1f%% (binomial P=%.2g) ; |effect| larger in tumour %.1f%% ; per-gene FDR<0.10: stronger in tumour %d, weaker %d",
      lab, 100 * bt$estimate, bt$p.value, 100 * mean(abs(G$mT) > abs(G$mN)), G[fdr < 0.10 & mT > mN, .N], G[fdr < 0.10 & mT < mN, .N])
  Pp <- merge(merge(Pp, pur, by = "pid", all.x = TRUE), epi, by = "pid", all.x = TRUE)
  msg("[%s] patients with tumour > normal %.0f%% (paired Wilcoxon P=%.2g) ; patient difference vs purity rho %+.3f (P=%.2g), vs epithelial rho %+.3f (P=%.2g)", lab,
      100 * Pp[, mean(mT > mN)], wilcox.test(Pp$mT, Pp$mN, paired = TRUE)$p.value, cor(Pp$d, Pp$purity, method = "s", use = "c"), cor.test(Pp$d, Pp$purity, method = "s", exact = FALSE)$p.value,
      cor(Pp$d, Pp$epiT, method = "s", use = "c"), cor.test(Pp$d, Pp$epiT, method = "s", exact = FALSE)$p.value)
  L <- rbind(A[, .(eg, pid, tissue = 1, y = aT)], A[, .(eg, pid, tissue = 0, y = aN)])
  fm <- lme4::lmer(y ~ tissue + (1 | eg) + (1 | pid) + (1 | eg:pid), data = L)
  cf <- summary(fm)$coefficients; b0 <- cf["(Intercept)", 1]; b1 <- cf["tissue", 1]; se <- cf["tissue", 2]
  msg("[%s] mixed model: normal mean %+.4f ; tumour - normal %+.4f [%+.4f, %+.4f] ; z = %.2f, P = %.2g ; tumour/normal %.3f", lab, b0, b1, b1 - 1.96 * se, b1 + 1.96 * se, b1 / se, 2 * pnorm(-abs(b1 / se)), (b0 + b1) / b0)
  fwrite(G[order(p)], file.path(O, sprintf("delta_eqtl_genes_%s.csv", lab))); fwrite(Pp, file.path(O, sprintf("delta_eqtl_patients_%s.csv", lab)))
  invisible(G) }
for (dd in c(FALSE, TRUE)) { sfx <- if (dd) "_dedup" else ""
  paired(QN, paste0("normal_map", sfx), dd); paired(QN[eg %chin% QT$eg], paste0("both_maps", sfx), dd); paired(QT, paste0("tumour_map", sfx), dd)
  if (nrow(QG)) { paired(QG, paste0("gtex_prostate", sfx), dd); paired(QG[eg %chin% QN$eg], paste0("gtex_and_normal", sfx), dd) } }
msg("DONE")
