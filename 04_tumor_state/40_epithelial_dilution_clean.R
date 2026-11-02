#!/usr/bin/env Rscript
# _clean (2026-10-06): input delta_ase_clean/paired_balanced_sites.csv.gz (56: overflow rows and wrap-risk genes removed), weights capped at 100 reads; outputs to delta_ase_clean / delta_eqtl_clean.
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

dc <- function(f) { d <- fread(f); d[, .(pid = sub("^[TN]", "", get(names(d)[1])), epi = Epithelial_cells)] }
eT <- dc(file.path(V, "decon_v2d/CPGEA_tumor_cell_fraction.tsv")); eN <- dc(file.path(V, "decon_v2d/CPGEA_normal_cell_fraction.tsv"))
rt <- function(f) { d <- fread(f); setnames(d, 1, "sid"); d[, pid := sub("^[TN]", "", sid)]; d }
tT <- rt(file.path(V, "decon/CPGEA_tumor_gene.tpm")); tN <- rt(file.path(V, "decon/CPGEA_normal_gene.tpm"))
mkg <- c("EPCAM", "KRT8", "KRT18", "CDH1", "CLDN4", "KRT19")
both <- rbind(tT[, c("pid", mkg), with = FALSE][, tis := "T"], tN[, c("pid", mkg), with = FALSE][, tis := "N"])
both[, mk := rowMeans(scale(log2(as.matrix(.SD) + 1))), .SDcols = mkg]
cv <- Reduce(function(a, b) merge(a, b, by = "pid"), list(eT[, .(pid, epiT = epi)], eN[, .(pid, epiN = epi)],
          both[tis == "T", .(pid, mkT = mk)], both[tis == "N", .(pid, mkN = mk)]))
cv[, `:=`(dEpi = epiT - epiN, dMk = mkT - mkN)]
msg("epithelial fraction median: tumour %.2f, normal %.2f ; delta median %.2f", median(cv$epiT), median(cv$epiN), median(cv$dEpi))

pat <- function(sub) sub[, .(sT = weighted.mean(sT, w), sN = weighted.mean(sN, w), n = .N), by = pid][, gain := sT - sN]
rat <- function(p) { bt <- replicate(2000, { q <- p[sample(.N, replace = TRUE)]; q[, weighted.mean(sT, n)] / q[, weighted.mean(sN, n)] })
  c(ratio = p[, weighted.mean(sT, n)] / p[, weighted.mean(sN, n)], lo = quantile(bt, .025, names = FALSE), hi = quantile(bt, .975, names = FALSE), n = nrow(p)) }

msg("\n(1) per-patient gain (sT - sN) regressed on delta epithelial content; intercept = gain at equal epithelial content")
for (b in c("all", "low", "mid", "high")) {
  p <- merge(pat(if (b == "all") A0 else A0[base == b]), cv, by = "pid")
  for (v in c("dEpi", "dMk")) { f <- summary(lm(as.formula(paste("gain ~", v)), data = p, weights = n))$coefficients
    rho <- cor.test(p$gain, p[[v]], method = "s", exact = FALSE)
    msg("  baseline %-4s %-4s: rho %+.3f (P=%.2g) ; intercept %+.4f (SE %.4f, P=%.2g) ; mean gain %+.4f ; slope %+.4f (P=%.2g)", b, v,
        rho$estimate, rho$p.value, f[1, 1], f[1, 2], f[1, 4], p[, weighted.mean(gain, n)], f[2, 1], f[2, 4]) } }

msg("\n(2) tumour/normal ratio in patient subsets"); SUB <- list()
sets <- list(all = cv$pid, `normal epi >= 0.6` = cv[epiN >= 0.6, pid], `normal epi >= 0.7` = cv[epiN >= 0.7, pid],
             `smallest |delta epi| tertile` = cv[abs(dEpi) <= quantile(abs(dEpi), 1/3), pid],
             `smallest |delta marker| tertile` = cv[abs(dMk) <= quantile(abs(dMk), 1/3), pid])
for (nm in names(sets)) for (b in c("all", "low", "high")) {
  sub <- A0[pid %chin% sets[[nm]]]; if (b != "all") sub <- sub[base == b]; r <- rat(pat(sub))
  SUB[[length(SUB) + 1]] <- data.table(subset = nm, baseline = b, ratio = r["ratio"], lo = r["lo"], hi = r["hi"], patients = r["n"]); msg("  %-32s baseline %-4s : %.3f [%.3f, %.3f] (patients %d)", nm, b, r["ratio"], r["lo"], r["hi"], r["n"]) }

sg <- fread(file.path(V, "decon_v2/signature_matrix_GSE185344v2.tsv")); setnames(sg, 1, "sym")
ct <- setdiff(names(sg), "sym"); sg[, spec := Epithelial_cells / rowSums(.SD), .SDcols = ct]
sp <- merge(unique(gtf[, .(eg, sym)], by = "eg"), sg[is.finite(spec), .(sym, spec)], by = "sym")[, .(spec = max(spec)), by = eg]
A1 <- merge(A0, sp, by = "eg")
msg("\n(3) epithelial specificity of host genes (epithelial share of summed cell-type means), median by baseline:")
print(unique(A1[, .(snp, base, spec)])[, .(snps = .N, median_spec = round(median(spec), 3)), by = base][order(base)])
A1[, stert := cut(spec, quantile(unique(A1[, .(snp, spec)])$spec, 0:3 / 3), include.lowest = TRUE, labels = c("low spec", "mid spec", "high spec"))]
msg("tumour/normal ratio by baseline x epithelial specificity:")
for (b in c("low", "mid", "high")) for (s in c("low spec", "mid spec", "high spec")) { r <- rat(pat(A1[base == b & stert == s]))
  msg("  baseline %-4s %-9s : %.3f [%.3f, %.3f]", b, s, r["ratio"], r["lo"], r["hi"]) }
fwrite(merge(pat(A0[base == "low"]), cv, by = "pid"), file.path(D, "dilution_lowbaseline_patients.csv")); msg("DONE")
fwrite(rbindlist(SUB), file.path(D, "subset_ratios_40.csv"))
