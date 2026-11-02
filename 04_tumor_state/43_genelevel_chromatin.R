#!/usr/bin/env Rscript
# 43: gene-level chromatin test (replaces the reporter-SNP position test of 42-A, since reporter SNPs are transcribed SNPs and not
# the regulatory variants): host-gene promoter (TSS +-2 kb) and +-50 kb neighbourhood overlap with TCGA PRAD ATAC peaks, and with
# PRAD-specific peaks (peak in <= 2 of 22 other cancer types). Set-up copied from 42:
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
gtf <- fread(cmd = "zcat /path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz | awk -F'\\t' '$3==\"gene\"' | cut -f1,4,5,7,9", header = FALSE, sep = "\t", col.names = c("gchr", "gs", "ge", "strand", "V1"))
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
rd <- function(t) fread(file.path(C, paste0(t, "_peakCalls.txt")), select = 1:3, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)][, .(chrom, s, e)]
P <- rd("PRAD"); setkey(P, chrom, s, e)
oth <- rbindlist(lapply(setdiff(types, "PRAD"), function(t) rd(t)[, type := t])); setkey(oth, chrom, s, e)
x <- foverlaps(P, oth, nomatch = 0L)[, .(n_other = uniqueN(type)), by = .(chrom, s = i.s, e = i.e)]
P <- merge(P, x, by = c("chrom", "s", "e"), all.x = TRUE)[is.na(n_other), n_other := 0L]
Ps <- P[n_other <= 2]; setkey(P, chrom, s, e); setkey(Ps, chrom, s, e)
msg("PRAD peaks %s ; PRAD-specific (<= 2 other types) %s", format(nrow(P), big.mark = ","), format(nrow(Ps), big.mark = ","))
gtf[, `:=`(eg = sub('.*gene_id "([^"]+)".*', "\\1", V1), chrom = as.character(gchr))][, tss := fifelse(strand == "+", gs, ge)]
G <- unique(gtf[eg %chin% unique(A0$eg), .(eg, chrom, tss, gs, ge)], by = "eg")
cnt <- function(Q, pk) { setkey(Q, chrom, s, e); foverlaps(Q, pk, nomatch = 0L)[, .N, by = eg] }
pr2 <- G[, .(eg, chrom, s = tss - 2000L, e = tss + 2000L)]; nb <- G[, .(eg, chrom, s = pmax(1L, gs - 50000L), e = ge + 50000L)]
G <- merge(G, cnt(copy(pr2), P)[, .(eg, prom_prad = N > 0)], by = "eg", all.x = TRUE)
G <- merge(G, cnt(copy(pr2), Ps)[, .(eg, prom_spec = N > 0)], by = "eg", all.x = TRUE)
G <- merge(G, cnt(copy(nb), Ps)[, .(eg, nb_spec = N)], by = "eg", all.x = TRUE)
G[is.na(prom_prad), prom_prad := FALSE][is.na(prom_spec), prom_spec := FALSE][is.na(nb_spec), nb_spec := 0L]
G[, nb_cls := fifelse(nb_spec == 0, "0", fifelse(nb_spec <= 2, "1-2", ">=3"))][, glen := ge - gs + 1]
msg("host genes %s ; promoter in PRAD peak %.1f%% ; promoter in PRAD-specific peak %.1f%% ; >= 1 PRAD-specific peak within +-50 kb %.1f%%",
    format(nrow(G), big.mark = ","), 100 * mean(G$prom_prad), 100 * mean(G$prom_spec), 100 * mean(G$nb_spec > 0))
A0 <- merge(A0, G[, .(eg, prom_prad, prom_spec, nb_cls)], by = "eg")
msg("\npatient-level tumour/normal ratio by host-gene chromatin class:")
for (b in c("low", "high")) {
  for (v in list(c("prom_prad", "promoter in PRAD peak"), c("prom_spec", "promoter in PRAD-specific peak"))) for (val in c(TRUE, FALSE)) {
    z <- A0[base == b & get(v[1]) == val]; msg("  baseline %-4s %-32s %-5s SNPs %6d : %s", b, v[2], val, uniqueN(z$snp), if (uniqueN(z$snp) >= 100) rat(z) else "too few") }
  for (k in c("0", "1-2", ">=3")) { z <- A0[base == b & nb_cls == k]; msg("  baseline %-4s PRAD-specific peaks +-50 kb = %-3s        SNPs %6d : %s", b, k, uniqueN(z$snp), if (uniqueN(z$snp) >= 100) rat(z) else "too few") } }
GG <- A0[base == "low", .(gain = weighted.mean(sT - sN, w), np = .N), by = eg]
GG <- merge(GG, G, by = "eg")
f <- summary(lm(gain ~ prom_spec + nb_spec + log(glen) + log(np), data = GG, weights = np))$coefficients
msg("\ngene-level (low-baseline SNPs) gain ~ promoter PRAD-specific + n PRAD-specific peaks +-50 kb + log length + log obs (%d genes):", nrow(GG))
print(round(f[-1, , drop = FALSE], 5))
fwrite(G, file.path(C, "hostgene_chromatin.csv")); msg("DONE")
