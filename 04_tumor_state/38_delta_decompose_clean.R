#!/usr/bin/env Rscript
# _clean (2026-10-06): input delta_ase_clean (56), weights capped at 100 reads.
# Where does the patient-level Delta-ASE increase (27-A, ratio 1.110) come from? Same filtered sites and leave-one-out orientation
# as 27-A; patient-level tumour/normal ratio (patient bootstrap) and tumour - normal difference within strata of:
#   (1) GTEx v8 prostate eGene status of the host gene: eGene (qval < 0.05) / tested non-eGene / not tested in GTEx;
#   (2) normal-tissue expression tertile of the host gene (log2 TPM, median over patients);
#   (3) baseline imbalance in normal tissue: |leave-one-out normal-only deviation| (other patients) tertile — distinguishes
#       amplification of imbalance already present in normal tissue from imbalance appearing in tumour.
suppressPackageStartupMessages(library(data.table)); setDTthreads(8); set.seed(38)
V <- "/path/to/data_root/ASE/code/v8_round2"; D <- file.path(V, "delta_ase_clean"); F <- file.path(V, "mapfilter")
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
B <- fread(file.path(D, "paired_balanced_sites.csv.gz"))
B[, `:=`(snp = paste(chrom, pos, ref, alt, sep = ":"), chrom = as.character(chrom), eg = sub("\\..*", "", eg))]
sd <- fread(file.path(F, "genomicSuperDups.txt.gz"), select = 2:4, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
bl <- fread(file.path(F, "hg38-blacklist.v2.bed.gz"), select = 1:3, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
ov <- function(reg) { q <- unique(B[, .(chrom, s = pos, e = pos)]); r <- unique(reg[, .(chrom, s = s + 1L, e)]); setkey(r, chrom, s, e)
  foverlaps(q, r, nomatch = 0L)[, unique(paste(chrom, i.s))] }
bad <- c(ov(sd), ov(bl))
gtf <- fread(cmd = "zcat /path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz | awk -F'\\t' '$3==\"gene\"' | cut -f9", header = FALSE, sep = "\t")
gtf[, `:=`(eg = sub('.*gene_id "([^"]+)".*', "\\1", V1), sym = sub('.*gene_name "([^"]+)".*', "\\1", V1))]
pseudo <- gtf[grepl('gene_biotype "[a-z_]*pseudogene"', V1), eg]
mono <- unique(B, by = c("pid", "snp"))[, .(devN = abs(sum(rnN) / sum(rnN + anN) - 0.5)), by = snp][devN > 0.35, snp]
B <- B[!(paste(chrom, pos) %chin% bad | (chrom == "6" & pos >= 28510120 & pos <= 33480577) | eg %chin% pseudo | snp %chin% mono)]
B <- unique(B, by = c("pid", "snp")); B[, `:=`(dT = rnT + anT, dN = rnN + anN)]
gt <- fread(file.path(V, "delta_eqtl/Prostate.v8.egenes.txt.gz"))[, .(eg = sub("\\..*", "", gene_id), qval)]
eN <- fread(file.path(V, "decon/CPGEA_normal_gene.tpm")); setnames(eN, 1, "sid")
ex <- data.table(sym = setdiff(names(eN), "sid"))[, med := sapply(sym, function(s) median(log2(eN[[s]] + 1)))]
ex <- merge(ex, unique(gtf[, .(eg, sym)], by = "sym"), by = "sym")[, .(med = max(med)), by = eg]
B <- merge(B, gt, by = "eg", all.x = TRUE); B <- merge(B, ex, by = "eg", all.x = TRUE)
B[, gtex := fifelse(is.na(qval), "not tested in GTEx", fifelse(qval < 0.05, "GTEx eGene", "GTEx tested, not eGene"))]
B[, npat := .N, by = snp]; A0 <- B[npat >= 6]
tot <- A0[, .(RT = sum(rnT), DT = sum(dT), RN = sum(rnN), DN = sum(dN)), by = snp]
A0 <- merge(A0, tot, by = "snp")
A0[, loo := 0.5 * ((RT - rnT) / (DT - dT) + (RN - rnN) / (DN - dN)) - 0.5][, looN := (RN - rnN) / (DN - dN) - 0.5]
A0 <- A0[abs(loo) >= 0.05]
A0[, `:=`(sN = sign(loo) * (rnN / dN - 0.5), sT = sign(loo) * (rnT / dT - 0.5), w = pmin(dN, dT, 100))]
A0[, base := cut(abs(looN), quantile(abs(looN), 0:3 / 3), include.lowest = TRUE, labels = c("low", "mid", "high"))]
A0[is.finite(med), etert := cut(med, quantile(med, 0:3 / 3, na.rm = TRUE), include.lowest = TRUE, labels = c("low expr", "mid expr", "high expr"))]
msg("LOO-oriented obs %s, SNPs %s, patients %d (should match 27-A: 1,612,697 / 168,676 / 117)", format(nrow(A0), big.mark = ","), format(uniqueN(A0$snp), big.mark = ","), uniqueN(A0$pid))
summ <- function(sub, lab) {
  p <- sub[, .(sT = weighted.mean(sT, w), sN = weighted.mean(sN, w), n = .N), by = pid]
  bt <- replicate(2000, { q <- p[sample(.N, replace = TRUE)]; q[, weighted.mean(sT, n)] / q[, weighted.mean(sN, n)] })
  data.table(stratum = lab, snps = uniqueN(sub$snp), obs = nrow(sub), normal = p[, weighted.mean(sN, n)], tumour = p[, weighted.mean(sT, n)],
             ratio = p[, weighted.mean(sT, n)] / p[, weighted.mean(sN, n)], lo = quantile(bt, .025), hi = quantile(bt, .975),
             diff = p[, weighted.mean(sT - sN, n)], pct_patients_T_gt_N = 100 * p[, mean(sT > sN)]) }
R <- rbind(summ(A0, "all"),
           rbindlist(lapply(c("GTEx eGene", "GTEx tested, not eGene", "not tested in GTEx"), function(g) summ(A0[gtex == g], paste("GTEx:", g)))),
           rbindlist(lapply(c("low expr", "mid expr", "high expr"), function(g) summ(A0[etert == g], paste("normal expression:", g)))),
           rbindlist(lapply(c("low", "mid", "high"), function(g) summ(A0[base == g], paste("baseline normal imbalance:", g)))),
           rbindlist(lapply(c("low expr", "mid expr", "high expr"), function(e) rbindlist(lapply(c("GTEx eGene", "GTEx tested, not eGene"), function(g)
             summ(A0[etert == e & gtex == g], sprintf("%s x %s", e, g)))))))
options(width = 200); print(R[, .(stratum, snps, normal = round(normal, 4), tumour = round(tumour, 4), ratio = round(ratio, 3), lo = round(lo, 3), hi = round(hi, 3), diff = round(diff, 4), pct_patients_T_gt_N = round(pct_patients_T_gt_N))])
fwrite(R, file.path(D, "delta_ase_decomposition.csv")); msg("DONE")
