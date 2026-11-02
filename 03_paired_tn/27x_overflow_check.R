#!/usr/bin/env Rscript
# 27x: impact check of 16-bit count overflow (negative or depth >= 32767; found by 37b) on H1 (A). Identical to 27 (A) except OVF=1 drops those rows right after loading.
suppressPackageStartupMessages({library(data.table); library(parallel)}); setDTthreads(4); set.seed(27)
O <- "/path/to/data_root/ASE/code/v8_round2/delta_ase"; F <- "/path/to/data_root/ASE/code/v8_round2/mapfilter"
NC <- as.integer(Sys.getenv("NC", "16")); NSWAP <- as.integer(Sys.getenv("NSWAP", "20"))
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

B <- fread(file.path(O, "paired_balanced_sites.csv.gz"))
B[, `:=`(snp = paste(chrom, pos, ref, alt, sep = ":"), chrom = as.character(chrom), eg = sub("\\..*", "", eg))]
OVF <- B[, rnT < 0 | anT < 0 | rnN < 0 | anN < 0 | rnT + anT >= 32767 | rnN + anN >= 32767]; msg("overflow rows %d (patients %d)", sum(OVF), uniqueN(B$pid[OVF])); OVROWS <- B[OVF, .(pid, snp)]
for (MODE in c("asis", "clean")) { B <- fread(file.path(O, "paired_balanced_sites.csv.gz")); B[, `:=`(snp = paste(chrom, pos, ref, alt, sep = ":"), chrom = as.character(chrom), eg = sub("\\..*", "", eg))]
if (MODE == "clean") B <- B[!(rnT < 0 | anT < 0 | rnN < 0 | anN < 0 | rnT + anT >= 32767 | rnN + anN >= 32767)]; msg("\n=== %s", MODE)
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
genes_of <- unique(B[, .(snp, eg)])
B <- unique(B[snp %chin% keepsnp], by = c("pid", "snp"))
B[, `:=`(dT = rnT + anT, dN = rnN + anN)]

A0 <- B[, npat := .N, by = snp][npat >= 6]
tot <- A0[, .(RT = sum(rnT), DT = sum(dT), RN = sum(rnN), DN = sum(dN)), by = snp]
A0 <- merge(A0, tot, by = "snp")[, loo := 0.5 * ((RT - rnT) / (DT - dT) + (RN - rnN) / (DN - dN)) - 0.5][abs(loo) >= 0.05]
A0[, `:=`(sN = sign(loo) * (rnN / dN - 0.5), sT = sign(loo) * (rnT / dT - 0.5))]
pat <- A0[, .(sN = weighted.mean(sN, pmin(dN, dT)), sT = weighted.mean(sT, pmin(dN, dT)), n = .N), by = pid]
msg("overflow obs left in A0: %d ; max weight %d", nrow(merge(A0[, .(pid, snp)], OVROWS, by = c("pid", "snp"))), A0[, max(pmin(dN, dT))])
bt <- replicate(2000, { p <- pat[sample(.N, replace = TRUE)]; p[, weighted.mean(sT, n)] / p[, weighted.mean(sN, n)] })
msg("(A) %s obs ; %d patients ; T/N %.4f [%.3f, %.3f] ; sT>sN %.0f%% ; P=%.2g", format(nrow(A0), big.mark=","), nrow(pat), pat[, weighted.mean(sT, n)] / pat[, weighted.mean(sN, n)], quantile(bt,.025), quantile(bt,.975), 100*pat[, mean(sT > sN)], wilcox.test(pat$sT, pat$sN, paired = TRUE)$p.value) }
