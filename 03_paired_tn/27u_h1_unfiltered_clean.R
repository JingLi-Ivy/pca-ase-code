#!/usr/bin/env Rscript
# 27u: H1 sensitivity on the cleaned data (56; depth weights capped at 100), with vs without the mappability filters of 27_clean.
# Lines below up to keepsnp are copied from 27_delta_ase_mapfilter_clean.R; summary (A) computed identically for both sets, same bootstrap seed.
# _clean (2026-10-06): input delta_ase_clean/paired_balanced_sites.csv.gz (56: overflow rows and wrap-risk genes removed), weights capped at 100 reads; outputs to delta_ase_clean / delta_eqtl_clean.
# Review follow-up for Delta-ASE: significant SNPs clustered in mapping-prone loci (LINC00910 17q21.31 duplication, HLA-A,
# pseudogenes), which WASP does not correct. Remove, before any testing:
#   (1) UCSC genomicSuperDups ; (2) ENCODE blacklist v2 ; (3) MHC chr6:28,510,120-33,480,577 ;
#   (4) SNPs whose host gene (any mapped gene) is a pseudogene (Ensembl 110 biotype) ;
#   (5) consistent near-monoallelic SNPs: pooled normal-tissue deviation |ref fraction - 0.5| > 0.35 across heterozygotes
#       (true imprinting is parent-of-origin, so its direction relative to ref/alt varies by patient; a consistent direction
#        across patients points to mapping or genotyping artefact).
# Then rerun: (A) patient-level LOO summary (as 19 part A); (B) per-SNP test (as 23a v2) with results collapsed to genes;
# (C) split-sample direction asymmetry (as 23c) with gene-clustered bootstrap CIs;
# (D) calibration: within-patient random swap of tumour/normal labels (one coin per patient, applied to all SNPs), NSWAP times,
#     rerunning (B) -> null number of significant SNPs/genes and null fraction tumour-stronger.
suppressPackageStartupMessages({library(data.table); library(parallel)}); setDTthreads(4); set.seed(27)
O <- "/path/to/data_root/ASE/code/v8_round2/delta_ase_clean"; F <- "/path/to/data_root/ASE/code/v8_round2/mapfilter"
NC <- as.integer(Sys.getenv("NC", "16")); NSWAP <- as.integer(Sys.getenv("NSWAP", "20"))
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

B <- fread(file.path(O, "paired_balanced_sites.csv.gz"))
B[, `:=`(snp = paste(chrom, pos, ref, alt, sep = ":"), chrom = as.character(chrom), eg = sub("\\..*", "", eg))]
sd <- fread(file.path(F, "genomicSuperDups.txt.gz"), select = 2:4, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
bl <- fread(file.path(F, "hg38-blacklist.v2.bed.gz"), select = 1:3, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
ov <- function(reg) { q <- unique(B[, .(chrom, s = pos, e = pos)]); r <- unique(reg[, .(chrom, s = s + 1L, e)]); setkey(r, chrom, s, e)
  unique(foverlaps(q, r, nomatch = 0L)[, .(chrom, pos = i.s)]) }
sdp <- ov(sd); blp <- ov(bl)
gtf <- fread(cmd = "zcat /path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz | awk -F'\\t' '$3==\"gene\"' | cut -f9", header = FALSE, sep = "\t")
gtf[, `:=`(eg = sub('.*gene_id "([^"]+)".*', "\\1", V1), bt = sub('.*gene_biotype "([^"]+)".*', "\\1", V1))]
pseudo <- gtf[grepl("pseudogene", bt), eg]
U <- unique(B, by = c("pid", "snp"))
mono <- U[, .(devN = abs(sum(rnN) / sum(rnN + anN) - 0.5)), by = snp][devN > 0.35, snp]
B[, `:=`(f_segdup = paste(chrom, pos) %chin% sdp[, paste(chrom, pos)], f_black = paste(chrom, pos) %chin% blp[, paste(chrom, pos)],
         f_mhc = chrom == "6" & pos >= 28510120 & pos <= 33480577, f_pseudo = eg %chin% pseudo, f_mono = snp %chin% mono)]
fl <- B[, .(f_segdup = any(f_segdup), f_black = any(f_black), f_mhc = any(f_mhc), f_pseudo = any(f_pseudo), f_mono = any(f_mono)), by = snp]
fl[, drop := f_segdup | f_black | f_mhc | f_pseudo | f_mono]
msg("SNPs %s ; removed: segdup %s, blacklist %s, MHC %s, pseudogene host %s, consistent near-monoallelic %s ; any %s (%.1f%%)",
    format(nrow(fl), big.mark = ","), format(sum(fl$f_segdup), big.mark = ","), format(sum(fl$f_black), big.mark = ","), format(sum(fl$f_mhc), big.mark = ","),
    format(sum(fl$f_pseudo), big.mark = ","), format(sum(fl$f_mono), big.mark = ","), format(sum(fl$drop), big.mark = ","), 100 * mean(fl$drop))
keepsnp <- fl[drop == FALSE, snp]
h1 <- function(X, lab) {
  X <- unique(X, by = c("pid", "snp"))[, `:=`(dT = rnT + anT, dN = rnN + anN)]
  A0 <- X[, npat := .N, by = snp][npat >= 6]
  tot <- A0[, .(RT = sum(rnT), DT = sum(dT), RN = sum(rnN), DN = sum(dN)), by = snp]
  A0 <- merge(A0, tot, by = "snp")[, loo := 0.5 * ((RT - rnT) / (DT - dT) + (RN - rnN) / (DN - dN)) - 0.5][abs(loo) >= 0.05]
  A0[, `:=`(sN = sign(loo) * (rnN / dN - 0.5), sT = sign(loo) * (rnT / dT - 0.5))]
  pat <- A0[, .(sN = weighted.mean(sN, pmin(dN, dT, 100)), sT = weighted.mean(sT, pmin(dN, dT, 100)), n = .N), by = pid]
  set.seed(1); bt <- replicate(2000, { p <- pat[sample(.N, replace = TRUE)]; p[, weighted.mean(sT, n)] / p[, weighted.mean(sN, n)] })
  msg("(%s) %s obs, %s SNPs, %d patients ; tumour/normal %.3f [%.3f, %.3f] ; patients with sT>sN %.0f%%, paired Wilcoxon P=%.2g", lab,
      format(nrow(A0), big.mark = ","), format(uniqueN(A0$snp), big.mark = ","), nrow(pat), pat[, weighted.mean(sT, n)] / pat[, weighted.mean(sN, n)],
      quantile(bt, .025), quantile(bt, .975), 100 * pat[, mean(sT > sN)], wilcox.test(pat$sT, pat$sN, paired = TRUE)$p.value) }
h1(copy(B), "unfiltered")
h1(B[snp %chin% keepsnp], "filtered")
msg("DONE")
