#!/usr/bin/env Rscript
# _clean (2026-10-06): input delta_ase_clean/paired_balanced_sites.csv.gz (56: overflow rows and wrap-risk genes removed), weights capped at 100 reads; outputs to delta_ase_clean / delta_eqtl_clean.
# 51: H4 (tumour gain tied to germline alleles) re-estimated over the same 20 random splits design as H3 (46),
# replacing the single split of 39 (2b). Data preparation identical to 39/46 (mappability-filtered paired balanced sites).
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
  H2 <- merge(S[half == 2], H1, by = "snp")
  H2[, `:=`(rT = rnT / dT - 0.5, rN = rnN / dN - 0.5, w = pmin(dN, dT, 100))]
  H2[, base := cut(abs(dN1), quantile(abs(dN1), 0:3 / 3), include.lowest = TRUE, labels = c("low", "mid", "high"))]
  rbindlist(lapply(c("dT1", "dN1"), function(col) rbindlist(lapply(c("all", "low", "high"), function(b) {
    z <- if (b == "all") H2 else H2[base == b]; o <- z[[col]]; z <- z[abs(o) >= 0.05]; o <- o[abs(o) >= 0.05]
    if (uniqueN(z$snp) < 200) return(NULL)
    z[, `:=`(sT = sign(o) * rT, sN = sign(o) * rN)]
    p <- z[, .(sT = weighted.mean(sT, w), sN = weighted.mean(sN, w), n = .N), by = pid]
    bt <- replicate(1000, { q <- p[sample(.N, replace = TRUE)]; q[, weighted.mean(sT, n)] / q[, weighted.mean(sN, n)] })
    data.table(seed = seed, orient = if (col == "dT1") "half-1 tumour" else "half-1 normal", base = b, snps = uniqueN(z$snp),
               ratio = p[, weighted.mean(sT, n)] / p[, weighted.mean(sN, n)], lo = quantile(bt, .025), hi = quantile(bt, .975), p_le1 = mean(bt <= 1)) })))) }
R <- rbindlist(parallel::mclapply(1:20, function(k) one_split(3900 + k), mc.cores = 10))
fwrite(R, file.path(D, "h4_splits.csv"))
msg("H4, 20 random splits (median across splits of per-split estimate / bootstrap CI / one-sided bootstrap P for ratio <= 1):")
print(R[, .(splits = .N, snps = median(snps), ratio = median(ratio), lo = median(lo), hi = median(hi), P = median(p_le1),
            ratio_min = min(ratio), ratio_max = max(ratio)), by = .(orient, base)])
msg("single split in 39 (seed 39): tumour-oriented all 1.170 [1.128, 1.229], weak 1.508, strong 1.122; normal-oriented all 1.074")
msg("DONE")
