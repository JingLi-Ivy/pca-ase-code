#!/usr/bin/env Rscript
# _clean (2026-10-06): input delta_ase_clean/paired_balanced_sites.csv.gz (56: overflow rows and wrap-risk genes removed), weights capped at 100 reads; outputs to delta_ase_clean / delta_eqtl_clean.
# 41: same question as 40 for WITHIN-epithelium state composition (basal vs luminal; tumour-derived vs benign-derived states),
# from the BayesPrism v2d state fractions. Set-up copied from 40:
# Is the "tumour-awakened" Delta-ASE gain (38/39: largest at sites with little imbalance in normal tissue) explained by epithelial
# dilution (normal epithelial fraction ~0.48 vs tumour ~0.78; epithelium-specific cis effects would look weak in stroma-rich normal)?
# Same filtered sites and leave-one-out orientation as 38; baseline = tertile of |LOO normal-only deviation|.
# (1) per-patient gain (weighted sT - sN) vs delta epithelial (tumour - normal; BayesPrism v2d, and an epithelial marker score):
#     regression intercept = gain expected at equal epithelial content in both tissues;
# (2) gain in patients whose normal tissue is epithelium-rich (>= 0.6, >= 0.7) and in the tertile with the smallest |delta epithelial|;
# (3) epithelial specificity of host genes (GSE185344 v2 signature: epithelial share of the summed cell-type means) by baseline group,
#     and the low-baseline gain by specificity tertile.
suppressPackageStartupMessages(library(data.table)); setDTthreads(8); set.seed(40)
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
B <- unique(B, by = c("pid", "snp")); B[, `:=`(dT = rnT + anT, dN = rnN + anN, pid = as.character(pid))]
B[, npat := .N, by = snp]; A0 <- B[npat >= 6]
tot <- A0[, .(RT = sum(rnT), DT = sum(dT), RN = sum(rnN), DN = sum(dN)), by = snp]
A0 <- merge(A0, tot, by = "snp")
A0[, loo := 0.5 * ((RT - rnT) / (DT - dT) + (RN - rnN) / (DN - dN)) - 0.5][, looN := (RN - rnN) / (DN - dN) - 0.5]
A0 <- A0[abs(loo) >= 0.05]
A0[, `:=`(sN = sign(loo) * (rnN / dN - 0.5), sT = sign(loo) * (rnT / dT - 0.5), w = pmin(dN, dT, 100))]
A0[, base := cut(abs(looN), quantile(abs(looN), 0:3 / 3), include.lowest = TRUE, labels = c("low", "mid", "high"))]

st <- fread(file.path(V, "decon_v2d/CPGEA_state_fraction.tsv"))
st[, `:=`(tis = substr(sample, 1, 1), pid = sub("^[TN]", "", sample))]
st[, `:=`(basal = Basal_benign + Basal_tumour, lum = Luminal_benign + Luminal_tumour + Proximal_Luminal_benign + Proximal_Luminal_tumour,
          tumlike = Basal_tumour + Luminal_tumour + Proximal_Luminal_tumour)]
st[, `:=`(epi = basal + lum)][, `:=`(basal_share = basal / epi, tum_share = tumlike / epi)]
cv <- merge(st[tis == "T", .(pid, bT = basal_share, tT = tum_share, eT = epi)], st[tis == "N", .(pid, bN = basal_share, tN = tum_share, eN = epi)], by = "pid")
cv[, `:=`(dBasal = bT - bN, dTum = tT - tN, dEpi = eT - eN)]
msg("median basal share of epithelium: tumour %.3f, normal %.3f ; tumour-derived share: tumour %.3f, normal %.3f",
    median(cv$bT), median(cv$bN), median(cv$tT), median(cv$tN))
pat <- function(sub) sub[, .(sT = weighted.mean(sT, w), sN = weighted.mean(sN, w), n = .N), by = pid][, gain := sT - sN]
rat <- function(p) { bt <- replicate(2000, { q <- p[sample(.N, replace = TRUE)]; q[, weighted.mean(sT, n)] / q[, weighted.mean(sN, n)] })
  c(ratio = p[, weighted.mean(sT, n)] / p[, weighted.mean(sN, n)], lo = quantile(bt, .025, names = FALSE), hi = quantile(bt, .975, names = FALSE), n = nrow(p)) }
msg("\n(1) per-patient gain ~ delta state composition (weighted); intercept = gain at equal composition")
for (b in c("low", "high")) { p <- merge(pat(A0[base == b]), cv, by = "pid")
  for (fm in c("gain ~ dBasal", "gain ~ dTum", "gain ~ dBasal + dTum + dEpi")) { f <- summary(lm(as.formula(fm), data = p, weights = n))$coefficients
    msg("  baseline %-4s %-28s intercept %+.4f (P=%.2g) ; mean gain %+.4f ; %s", b, fm, f[1, 1], f[1, 4], p[, weighted.mean(gain, n)],
        paste(sprintf("%s %+.4f (P=%.2g)", rownames(f)[-1], f[-1, 1], f[-1, 4]), collapse = " ; ")) } }
msg("\n(2) tumour/normal ratio in patient subsets"); SUB <- list()
sets <- list(all = cv$pid, `lowest normal basal share (tertile)` = cv[bN <= quantile(bN, 1/3), pid],
             `smallest |delta basal share| tertile` = cv[abs(dBasal) <= quantile(abs(dBasal), 1/3), pid],
             `smallest |delta tumour-derived share| tertile` = cv[abs(dTum) <= quantile(abs(dTum), 1/3), pid])
for (nm in names(sets)) for (b in c("low", "high")) { r <- rat(pat(A0[pid %chin% sets[[nm]] & base == b]))
  SUB[[length(SUB) + 1]] <- data.table(subset = nm, baseline = b, ratio = r["ratio"], lo = r["lo"], hi = r["hi"], patients = r["n"]); msg("  %-46s baseline %-4s : %.3f [%.3f, %.3f] (patients %d)", nm, b, r["ratio"], r["lo"], r["hi"], r["n"]) }
res <- readRDS(file.path(V, "decon_v2d/bayesprism_v2c_result.rds"))
phi <- tryCatch(res@prism@phi_cellState@phi, error = function(e) NULL)
if (!is.null(phi)) {
  L <- colMeans(phi[grep("^(Luminal|Proximal_Luminal)", rownames(phi)), , drop = FALSE]); Bs <- colMeans(phi[grep("^Basal", rownames(phi)), , drop = FALSE])
  ls <- data.table(sym = names(L), lspec = L / (L + Bs))[is.finite(lspec)]
  ls <- merge(unique(gtf[, .(eg, sym)], by = "eg"), ls, by = "sym")[, .(lspec = max(lspec)), by = eg]
  A1 <- merge(A0, ls, by = "eg"); msg("\n(3) luminal specificity L/(L+B), host genes in the BayesPrism marker set: %s SNPs", format(uniqueN(A1$snp), big.mark = ","))
  print(unique(A1[, .(snp, base, lspec)])[, .(snps = .N, median_luminal_spec = round(median(lspec), 3)), by = base][order(base)])
  A1[, lt := cut(lspec, quantile(unique(A1[, .(snp, lspec)])$lspec, 0:3 / 3), include.lowest = TRUE, labels = c("basal-leaning", "intermediate", "luminal-leaning"))]
  for (b in c("low", "high")) for (s in c("basal-leaning", "intermediate", "luminal-leaning")) { r <- rat(pat(A1[base == b & lt == s]))
    msg("  baseline %-4s %-16s : %.3f [%.3f, %.3f]", b, s, r["ratio"], r["lo"], r["hi"]) }
} else msg("(3) state reference profiles not found in the BayesPrism object")
msg("DONE")
fwrite(rbindlist(SUB), file.path(D, "subset_ratios_41.csv"))
