#!/usr/bin/env Rscript
# Mechanism leads for "tumour-awakened" germline cis effects (38-41: Delta-ASE gain concentrated at sites with weak imbalance in
# normal tissue). Same filtered sites, leave-one-out orientation and baseline tertiles as 38/40.
# (A) chromatin: TCGA ATAC-seq cancer-type peak sets (Corces et al. 2018 Science; hg38, 501 bp). A SNP is "in PRAD peak" if inside a
#     PRAD peak; PRAD-specific if the position is a peak in <= 2 of the other 22 cancer types. No normal-prostate ATAC is available
#     here, so this asks whether awakened sites are open in prostate cancer (and prostate-specifically), not whether they newly opened.
# (B) TF: ADASTRA v6.1 allele-specific binding at the SNP (any TF, FDR < 0.05; allele-matched). Per TF with >= 50 ASB SNPs: gain at
#     SNPs where the TF-preferred allele is the oriented high allele; across TFs, gain vs the TF's tumour/normal expression change.
suppressPackageStartupMessages(library(data.table)); setDTthreads(8); set.seed(42)
V <- "/path/to/data_root/ASE/code/v8_round2"; D <- file.path(V, "delta_ase"); F <- file.path(V, "mapfilter"); C <- file.path(V, "chromatin")
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
B <- fread(file.path(D, "paired_balanced_sites.csv.gz"))
B[, `:=`(snp = paste(chrom, pos, ref, alt, sep = ":"), chrom = as.character(chrom), eg = sub("\\..*", "", eg))]
sd <- fread(file.path(F, "genomicSuperDups.txt.gz"), select = 2:4, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
bl <- fread(file.path(F, "hg38-blacklist.v2.bed.gz"), select = 1:3, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
ov <- function(reg, q = unique(B[, .(chrom, s = pos, e = pos)])) { r <- unique(reg[, .(chrom, s = s + 1L, e)]); setkey(r, chrom, s, e)
  foverlaps(q, r, nomatch = 0L)[, unique(paste(chrom, i.s))] }
bad <- c(ov(sd), ov(bl))
gtf <- fread(cmd = "zcat /path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz | awk -F'\\t' '$3==\"gene\"' | cut -f9", header = FALSE, sep = "\t")
pseudo <- gtf[grepl('gene_biotype "[a-z_]*pseudogene"', V1), sub('.*gene_id "([^"]+)".*', "\\1", V1)]
mono <- unique(B, by = c("pid", "snp"))[, .(devN = abs(sum(rnN) / sum(rnN + anN) - 0.5)), by = snp][devN > 0.35, snp]
B <- B[!(paste(chrom, pos) %chin% bad | (chrom == "6" & pos >= 28510120 & pos <= 33480577) | eg %chin% pseudo | snp %chin% mono)]
B <- unique(B, by = c("pid", "snp")); B[, `:=`(dT = rnT + anT, dN = rnN + anN)]
B[, npat := .N, by = snp]; A0 <- B[npat >= 6]
tot <- A0[, .(RT = sum(rnT), DT = sum(dT), RN = sum(rnN), DN = sum(dN)), by = snp]
A0 <- merge(A0, tot, by = "snp")
A0[, loo := 0.5 * ((RT - rnT) / (DT - dT) + (RN - rnN) / (DN - dN)) - 0.5][, looN := (RN - rnN) / (DN - dN) - 0.5]
A0 <- A0[abs(loo) >= 0.05]
A0[, `:=`(sN = sign(loo) * (rnN / dN - 0.5), sT = sign(loo) * (rnT / dT - 0.5), w = pmin(dN, dT))]
A0[, base := cut(abs(looN), quantile(abs(looN), 0:3 / 3), include.lowest = TRUE, labels = c("low", "mid", "high"))]
S <- A0[, .(gain = weighted.mean(sT - sN, w), np = .N, dep = sum(dT + dN), base = base[1], chrom = chrom[1], pos = pos[1],
            hi_ref = (RT[1] + RN[1]) / (DT[1] + DN[1]) > 0.5), by = snp]
pat <- function(sub) sub[, .(sT = weighted.mean(sT, w), sN = weighted.mean(sN, w), n = .N), by = pid]
rat <- function(sub) { p <- pat(sub); bt <- replicate(1000, { q <- p[sample(.N, replace = TRUE)]; q[, weighted.mean(sT, n)] / q[, weighted.mean(sN, n)] })
  sprintf("%.3f [%.3f, %.3f]", p[, weighted.mean(sT, n)] / p[, weighted.mean(sN, n)], quantile(bt, .025), quantile(bt, .975)) }

types <- sub("_peakCalls.txt", "", list.files(C, pattern = "_peakCalls.txt$"))
q <- S[, .(chrom, s = pos, e = pos)]
hit <- sapply(types, function(t) { pk <- fread(file.path(C, paste0(t, "_peakCalls.txt")), select = 1:3, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
  paste(S$chrom, S$pos) %chin% ov(pk, q) })
S[, `:=`(prad = hit[, "PRAD"], n_other = rowSums(hit[, setdiff(types, "PRAD"), drop = FALSE]))]
S[, atac := fifelse(!prad, "not in PRAD peak", fifelse(n_other <= 2, "PRAD-specific peak", "PRAD peak shared"))]
msg("(A) SNPs %s ; in PRAD ATAC peak %.1f%% (PRAD-specific %.2f%%)", format(nrow(S), big.mark = ","), 100 * mean(S$prad), 100 * mean(S$atac == "PRAD-specific peak"))
print(S[, .(snps = .N, pct_in_PRAD_peak = round(100 * mean(prad), 2), pct_PRAD_specific = round(100 * mean(atac == "PRAD-specific peak"), 3)), by = base][order(base)])
A0 <- merge(A0, S[, .(snp, atac, prad)], by = "snp")
msg("  patient-level tumour/normal ratio, low-baseline SNPs:")
for (a in c("not in PRAD peak", "PRAD peak shared", "PRAD-specific peak")) { z <- A0[base == "low" & atac == a]
  msg("    %-20s SNPs %6d : %s", a, uniqueN(z$snp), if (uniqueN(z$snp) >= 100) rat(z) else "too few") }
msg("  patient-level ratio, high-baseline SNPs in / not in PRAD peak: %s / %s", rat(A0[base == "high" & prad == TRUE]), rat(A0[base == "high" & prad == FALSE]))
f <- summary(lm(gain ~ prad + base + log(dep), data = S, weights = np))$coefficients
msg("  SNP-level gain ~ PRAD peak + baseline + log depth: PRAD peak %+.5f (P=%.2g)", f["pradTRUE", 1], f["pradTRUE", 4])
fl <- summary(lm(gain ~ prad + log(dep), data = S[base == "low"], weights = np))$coefficients
msg("  within low baseline: PRAD peak %+.5f (P=%.2g) ; mean low-baseline gain %+.5f", fl["pradTRUE", 1], fl["pradTRUE", 4], S[base == "low", weighted.mean(gain, np)])

A <- fread(file.path(V, "cistrans_validate/adastra_paired_sites.tsv"), header = FALSE, col.names = c("mn", "chr", "pos", "ref", "alt", "fdr_ref", "fdr_alt"))
tm <- fread("/path/to/data_root/ASE/code/adastra_tf_map.csv"); A <- merge(A, tm[, .(mn, tf = sym)], by = "mn")
A[, `:=`(snp = paste(sub("^chr", "", chr), pos, ref, alt, sep = ":"), pref_ref = fifelse(is.na(fdr_alt), TRUE, fifelse(is.na(fdr_ref), FALSE, fdr_ref < fdr_alt)))]
A <- unique(A[snp %chin% S$snp], by = c("snp", "tf"))
S[, asb := snp %chin% A$snp]
msg("\n(B) SNPs with any ASB: %d (%.2f%%) ; by baseline:", S[asb == TRUE, .N], 100 * mean(S$asb))
print(S[, .(snps = .N, pct_asb = round(100 * mean(asb), 2), gain_asb = round(weighted.mean(gain[asb], np[asb]), 5), gain_no = round(weighted.mean(gain[!asb], np[!asb]), 5)), by = base][order(base)])
A <- merge(A, S[, .(snp, gain, np, base, hi_ref)], by = "snp")[, concord := pref_ref == hi_ref]
rt <- function(f) { d <- fread(f); setnames(d, 1, "sid"); d[, pid := sub("^[TN]", "", sid)]; d }
eT <- rt(file.path(V, "decon/CPGEA_tumor_gene.tpm")); eN <- rt(file.path(V, "decon/CPGEA_normal_gene.tpm")); pp <- intersect(eT$pid, eN$pid)
tfs <- intersect(intersect(unique(A$tf), names(eT)), names(eN))
lfc <- data.table(tf = tfs, lfc = colMeans(log2(as.matrix(eT[match(pp, pid), ..tfs]) + 1) - log2(as.matrix(eN[match(pp, pid), ..tfs]) + 1)))
low_ref <- S[base == "low", weighted.mean(gain, np)]
TF <- A[, .(n = .N, n_conc = sum(concord), gain_conc = if (sum(concord) >= 5) weighted.mean(gain[concord], np[concord]) else NA_real_,
            gain_disc = if (sum(!concord) >= 5) weighted.mean(gain[!concord], np[!concord]) else NA_real_,
            p = tryCatch(wilcox.test(gain[concord], gain[!concord])$p.value, error = function(e) NA_real_)), by = tf][n >= 50]
TF <- merge(TF, lfc, by = "tf", all.x = TRUE)[, fdr := p.adjust(p, "BH")]
msg("TFs with >= 50 ASB SNPs: %d ; concordant vs discordant gain FDR<0.10: %d", nrow(TF), TF[fdr < 0.10, .N])
print(TF[order(p)][1:10, .(tf, n, n_conc, gain_conc = round(gain_conc, 4), gain_disc = round(gain_disc, 4), p = signif(p, 2), fdr = signif(fdr, 2), lfc = round(lfc, 2))])
ct <- cor.test(TF$gain_conc - TF$gain_disc, TF$lfc, method = "s", exact = FALSE)
msg("across TFs: (concordant - discordant gain) vs TF tumour/normal log2FC: rho %+.3f (P=%.2g, %d TFs)", ct$estimate, ct$p.value, TF[is.finite(lfc) & is.finite(gain_conc - gain_disc), .N])
fwrite(S, file.path(C, "awakened_snp_annotation.csv.gz")); fwrite(TF, file.path(C, "awakened_tf_asb.csv")); msg("DONE")
