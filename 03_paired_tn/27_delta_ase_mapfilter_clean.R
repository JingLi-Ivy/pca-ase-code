#!/usr/bin/env Rscript
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
R1 <- fread(file.path(O, "delta_ase_per_snp_v2.csv.gz"))[sig == TRUE]
msg("of the 63 v2 significant SNPs, removed by the filters: %d", R1[snp %in% fl[drop == TRUE, snp], .N])
keepsnp <- fl[drop == FALSE, snp]
genes_of <- unique(B[, .(snp, eg)])
B <- unique(B[snp %chin% keepsnp], by = c("pid", "snp"))
B[, `:=`(dT = rnT + anT, dN = rnN + anN)]

A0 <- B[, npat := .N, by = snp][npat >= 6]
tot <- A0[, .(RT = sum(rnT), DT = sum(dT), RN = sum(rnN), DN = sum(dN)), by = snp]
A0 <- merge(A0, tot, by = "snp")[, loo := 0.5 * ((RT - rnT) / (DT - dT) + (RN - rnN) / (DN - dN)) - 0.5][abs(loo) >= 0.05]
A0[, `:=`(sN = sign(loo) * (rnN / dN - 0.5), sT = sign(loo) * (rnT / dT - 0.5))]
pat <- A0[, .(sN = weighted.mean(sN, pmin(dN, dT, 100)), sT = weighted.mean(sT, pmin(dN, dT, 100)), n = .N, purity = purity[1]), by = pid]
bt <- replicate(2000, { p <- pat[sample(.N, replace = TRUE)]; p[, weighted.mean(sT, n)] / p[, weighted.mean(sN, n)] })
msg("\n(A) filtered: %s obs, %s SNPs, %d patients ; tumour/normal %.3f [%.3f, %.3f] ; patients with sT>sN %.0f%%, paired Wilcoxon P=%.2g",
    format(nrow(A0), big.mark = ","), format(uniqueN(A0$snp), big.mark = ","), nrow(pat), pat[, weighted.mean(sT, n)] / pat[, weighted.mean(sN, n)],
    quantile(bt, .025), quantile(bt, .975), 100 * pat[, mean(sT > sN)], wilcox.test(pat$sT, pat$sN, paired = TRUE)$p.value)
fwrite(pat, file.path(O, "mapfilter_A_per_patient.csv"))

mono2 <- function(r, d) r == 0 | r == d
Bi <- B[(mono2(rnT, dT) & mono2(rnN, dN)) == FALSE][, npat := .N, by = snp][npat >= 6]
fit1 <- function(z) {
  f <- tryCatch(suppressWarnings(glm(cbind(h, d - h) ~ factor(pid) + tissue, family = quasibinomial, data = z)), error = function(e) NULL)
  if (is.null(f) || f$df.residual < 1) return(list(b = NA_real_, p = NA_real_))
  phi <- sum(residuals(f, "pearson")^2) / f$df.residual; s <- summary(f, dispersion = max(1, phi))$coef
  if (!"tissue" %in% rownames(s)) return(list(b = NA_real_, p = NA_real_))
  list(b = s["tissue", 1], p = 2 * pt(-abs(s["tissue", 3]), f$df.residual)) }
pertest <- function(D, swap = NULL) {
  D <- copy(D)
  if (!is.null(swap)) { sw <- D$pid %in% swap; tmp <- D$rnT[sw]; D$rnT[sw] <- D$rnN[sw]; D$rnN[sw] <- tmp
    tmp <- D$anT[sw]; D$anT[sw] <- D$anN[sw]; D$anN[sw] <- tmp; D[, `:=`(dT = rnT + anT, dN = rnN + anN)] }
  S <- D[, .(avg = 0.5 * (sum(rnT) / sum(dT) + sum(rnN) / sum(dN)) - 0.5), by = snp][abs(avg) >= 0.05]
  X <- merge(D[snp %chin% S$snp], S[, .(snp, dir = sign(avg))], by = "snp")
  X[, `:=`(hT = ifelse(dir > 0, rnT, anT), hN = ifelse(dir > 0, rnN, anN))]
  L <- rbind(X[, .(snp, pid, tissue = 1L, h = hT, d = dT)], X[, .(snp, pid, tissue = 0L, h = hN, d = dN)]); setkey(L, snp)
  ch <- split(unique(L$snp), cut(seq_along(unique(L$snp)), NC, labels = FALSE))
  R <- rbindlist(mclapply(ch, function(ss) L[.(ss), fit1(.SD), by = snp], mc.cores = NC))[is.finite(p)]
  R <- merge(R, X[, .(fT = sum(hT) / sum(dT) - 0.5, fN = sum(hN) / sum(dN) - 0.5), by = snp], by = "snp")
  R[, fdr := p.adjust(p, "BH")][, sig := fdr < 0.05 & abs(fT - fN) >= 0.05]
  R[, cls := fcase(sig & b > 0, "tumor_stronger", sig & b < 0, "normal_stronger", default = "other")] }
genecount <- function(R) { g <- merge(R[sig == TRUE, .(snp, cls)], genes_of, by = "snp", allow.cartesian = TRUE)
  g <- g[, .(cls = if (all(cls == "tumor_stronger")) "tumor_stronger" else if (all(cls == "normal_stronger")) "normal_stronger" else "mixed"), by = eg]
  c(snp_T = R[cls == "tumor_stronger", .N], snp_N = R[cls == "normal_stronger", .N], gene_T = g[cls == "tumor_stronger", .N], gene_N = g[cls == "normal_stronger", .N],
    frac_T_p001 = R[p < 0.001, mean(b > 0)], n_p001 = R[p < 0.001, .N]) }
Rf <- pertest(Bi); ob <- genecount(Rf)
msg("\n(B) filtered per-SNP: tested %s ; significant SNPs tumour-stronger %d / normal-stronger %d ; genes %d / %d ; P<0.001 tests %d, tumour-stronger %.3f",
    format(nrow(Rf), big.mark = ","), ob["snp_T"], ob["snp_N"], ob["gene_T"], ob["gene_N"], ob["n_p001"], ob["frac_T_p001"])
sg <- merge(Rf[sig == TRUE], genes_of, by = "snp")[, .(snps = uniqueN(snp), cls = paste(unique(cls), collapse = "/"), med_fN = median(fN), med_fT = median(fT)), by = eg]
fwrite(Rf[order(p)], file.path(O, "mapfilter_B_per_snp.csv.gz")); fwrite(sg, file.path(O, "mapfilter_B_sig_genes.csv"))

pids <- unique(Bi$pid)
NUL <- rbindlist(lapply(seq_len(NSWAP), function(k) { set.seed(900 + k); sw <- pids[rbinom(length(pids), 1, 0.5) == 1]
  x <- genecount(pertest(Bi, swap = sw)); msg("  swap %d: SNPs T/N %d/%d, genes T/N %d/%d, frac tumour-stronger at P<0.001 %.3f (n=%d)",
    k, x["snp_T"], x["snp_N"], x["gene_T"], x["gene_N"], x["frac_T_p001"], x["n_p001"]); as.data.table(as.list(x)) }))
fwrite(NUL, file.path(O, "mapfilter_D_swap_null.csv"))
msg("(D) label-swap null (%d swaps): significant SNPs median %.0f (max %.0f) vs observed %d ; genes median %.0f (max %.0f) vs observed %d ; frac tumour-stronger at P<0.001: null %.3f +- %.3f vs observed %.3f",
    NSWAP, median(NUL$snp_T + NUL$snp_N), max(NUL$snp_T + NUL$snp_N), ob["snp_T"] + ob["snp_N"], median(NUL$gene_T + NUL$gene_N), max(NUL$gene_T + NUL$gene_N),
    ob["gene_T"] + ob["gene_N"], mean(NUL$frac_T_p001), sd(NUL$frac_T_p001), ob["frac_T_p001"])

Bs <- Bi[npat >= 8][, half := sample(rep(1:2, length.out = .N)), by = snp]
run <- function(oh) {
  S <- Bs[half == oh, .(avg = 0.5 * (sum(rnT) / sum(dT) + sum(rnN) / sum(dN)) - 0.5), by = snp][abs(avg) >= 0.05]
  X <- merge(Bs[half != oh & snp %chin% S$snp], S[, .(snp, dir = sign(avg))], by = "snp")
  X[, `:=`(hT = ifelse(dir > 0, rnT, anT), hN = ifelse(dir > 0, rnN, anN))]; X <- X[, if (.N >= 3) .SD, by = snp]
  L <- rbind(X[, .(snp, pid, tissue = 1L, h = hT, d = dT)], X[, .(snp, pid, tissue = 0L, h = hN, d = dN)]); setkey(L, snp)
  ch <- split(unique(L$snp), cut(seq_along(unique(L$snp)), NC, labels = FALSE))
  R <- rbindlist(mclapply(ch, function(ss) L[.(ss), fit1(.SD), by = snp], mc.cores = NC))[is.finite(p)]
  merge(R, X[, .(fT = sum(hT) / sum(dT) - 0.5, fN = sum(hN) / sum(dN) - 0.5), by = snp], by = "snp") }
SP <- rbind(run(1L), run(2L)); SP <- merge(SP, unique(Bs[, .(snp, eg)], by = "snp"), by = "snp")
msg("\n(C) filtered split-sample: %s held-out tests", format(nrow(SP), big.mark = ","))
for (th in c(1, 0.1, 0.01, 0.001)) { z <- SP[p < th]; if (nrow(z) < 10) next
  gs <- split(seq_len(nrow(z)), z$eg); bb <- replicate(2000, { i <- unlist(gs[sample(length(gs), replace = TRUE)], use.names = FALSE); mean(z$b[i] > 0) })
  msg("  P<%-6g n=%6d genes=%5d  tumour-stronger %.3f  gene-bootstrap [%.3f, %.3f]", th, nrow(z), length(gs), mean(z$b > 0), quantile(bb, .025), quantile(bb, .975)) }
gs <- split(seq_len(nrow(SP)), SP$eg); rb <- replicate(2000, { i <- unlist(gs[sample(length(gs), replace = TRUE)], use.names = FALSE); mean(SP$fT[i]) / mean(SP$fN[i]) })
msg("  held-out tumour/normal %.3f, gene-bootstrap [%.3f, %.3f]", SP[, mean(fT) / mean(fN)], quantile(rb, .025), quantile(rb, .975))
fwrite(SP, file.path(O, "mapfilter_C_split.csv.gz")); msg("DONE")
