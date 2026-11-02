#!/usr/bin/env Rscript
# 67: regenerate Supplemental Tables S14 and S17 for v10 (2026-10-10) from result files; no hand-copied numbers.
#  S14 = per-SNP paired tumour-normal allelic test (27_delta_ase_mapfilter_clean, part B): quasi-binomial GLM per SNP,
#        reads on the population-preferred allele ~ patient (fixed) + tissue, dispersion floored at 1; b = tissue log-odds (tumour vs normal).
#  S17 = transcriptome-wide cis-by-trans, slope statistic (26e), with covariate models, 1,000 whole-patient permutations,
#        measurement-precision bias (49, 50), prostate TF sets (52) and interaction-eQTL specificity (28f).
#  The previous S17 (correlation-difference statistic; ADASTRA MH-OR from an rsID merge error) is superseded.
suppressPackageStartupMessages(library(data.table))
V <- "/path/to/data_root/ASE/code/v8_round2"; O <- file.path(V, "supp_data"); setwd(V)

R <- fread("delta_ase_clean/mapfilter_B_per_snp.csv.gz")
S <- fread("delta_ase_clean/paired_balanced_sites.csv.gz", select = c("pid", "eg", "chrom", "pos", "ref", "alt"))
S[, snp := paste(chrom, pos, ref, alt, sep = ":")][, eg := sub("\\..*", "", eg)]
S <- S[snp %chin% R$snp]
np <- unique(S[, .(snp, pid)])[, .(n_patients = .N), by = snp]
ge <- unique(S[, .(snp, eg)])
gtf <- fread(cmd = "zcat /path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz | awk -F'\\t' '$3==\"gene\"' | cut -f9", header = FALSE, sep = "\t")
gtf[, `:=`(eg = sub('.*gene_id "([^"]+)".*', "\\1", V1), sym = fifelse(grepl("gene_name", V1), sub('.*gene_name "([^"]+)".*', "\\1", V1), NA_character_))]
ge <- merge(ge, unique(gtf[, .(eg, sym)]), by = "eg", all.x = TRUE)
ge <- ge[, .(ensembl_gene = paste(sort(unique(eg)), collapse = ";"), gene_symbol = paste(sort(unique(na.omit(sym))), collapse = ";")), by = snp]
T14 <- merge(merge(R, np, by = "snp", all.x = TRUE), ge, by = "snp", all.x = TRUE)
T14 <- T14[order(p), .(SNP = snp, ensembl_gene, gene_symbol, n_patients, tumor_deviation = round(fT, 5), normal_deviation = round(fN, 5),
                       difference = round(fT - fN, 5), tissue_log_odds = signif(b, 5), P = signif(p, 4), FDR = signif(fdr, 4), significant = sig,
                       class = fifelse(cls == "other", "", cls))]
stopifnot(nrow(T14) == 101053, T14[significant == TRUE, .N] == 23, T14[class == "tumor_stronger", .N] == 19, !anyNA(T14$n_patients))
fwrite(T14, file.path(O, "TableS14_per_SNP_paired_tumor_normal_test.csv"))
cat(sprintf("S14: %d SNPs, %d significant (%d tumor-stronger, %d normal-stronger); patients per SNP median %d\n",
            nrow(T14), T14[significant == TRUE, .N], T14[class == "tumor_stronger", .N], T14[class == "normal_stronger", .N], as.integer(median(T14$n_patients))))

rows <- list(); add <- function(sec, met, val, lo = NA, hi = NA, p = NA, n = NA, src = "") rows[[length(rows) + 1]] <<- data.table(Section = sec, Metric = met, Value = val, CI_lower = lo, CI_upper = hi, P = p, N = n, Source = src)
f4 <- function(x) sprintf("%+.4f", x)
CM <- fread("cistrans_slope_26e/tn_slope_by_covariate_model.csv"); NL <- fread("cistrans_slope_26e/tn_slope_strict_null.csv")
lab <- c(`ai ~ ld` = "Depth only (primary)", `ai ~ ld + epi` = "+ epithelial marker score", `ai ~ ld + pro` = "+ proliferation score",
         `ai ~ ld + epi + pro` = "+ epithelial marker score + proliferation score", `ai ~ ld + bp_epi` = "+ BayesPrism epithelial fraction")
pm <- c(`ai ~ ld` = "r0", `ai ~ ld + bp_epi` = "r_bp")
for (i in seq_len(nrow(CM))) { m <- CM$model[i]; sec <- sprintf("A. Mean slope-difference z by covariate model: %s", lab[[m]])
  nn <- sprintf("%d genes; %s gene x factor tests", CM$genes[i], format(CM$tests[i], big.mark = ","))
  add(sec, "Tumor", f4(CM$tumor[i]), f4(CM$tumor_lo[i]), f4(CM$tumor_hi[i]), NA, nn, "26e (gene bootstrap)")
  add(sec, "Normal", f4(CM$normal[i]), f4(CM$normal_lo[i]), f4(CM$normal_hi[i]), NA, nn, "26e (gene bootstrap)")
  add(sec, "Tumor minus normal", f4(CM$delta[i]), f4(CM$delta_lo[i]), f4(CM$delta_hi[i]), signif(CM$delta_P[i], 2), nn, "26e (gene bootstrap P)")
  if (m %in% names(pm)) { z <- NL[model == pm[[m]]]
    for (s in c("tumor", "normal", "delta")) { ob <- CM[[s]][i]; nu <- z[[s]]
      add(sec, sprintf("%s: whole-patient permutation (null mean +- SD)", c(tumor = "Tumor", normal = "Normal", delta = "Tumor minus normal")[[s]]),
          sprintf("%s +- %.4f", f4(mean(nu)), sd(nu)), NA, NA, signif((sum(abs(nu - mean(nu)) >= abs(ob - mean(nu))) + 1) / (length(nu) + 1), 2), sprintf("%d permutations", length(nu)), "26e strict null") } } }
P1 <- fread("reliability_bias/het_hom_snp_depth.csv.gz")
for (tis in c("tumour", "normal")) { q <- P1[tissue == tis]
  add("B. Measurement precision of |a| by genotype group (per-gene medians)", sprintf("%s: reporter SNPs per patient-gene, heterozygotes vs homozygotes", tis),
      sprintf("%.1f vs %.1f", q[het == TRUE, median(nsnp)], q[het == FALSE, median(nsnp)]), NA, NA, NA, sprintf("%d genes", uniqueN(q$gene)), "49")
  add("B. Measurement precision of |a| by genotype group (per-gene medians)", sprintf("%s: aggregate depth, heterozygotes vs homozygotes", tis),
      sprintf("%.0f vs %.0f", q[het == TRUE, median(dep)], q[het == FALSE, median(dep)]), NA, NA, NA, sprintf("%d genes", uniqueN(q$gene)), "49") }
SR <- fread("reliability_bias/simulation_reps.csv")[, .(fz = mean(fz), fz_lo = quantile(fz, .025), fz_hi = quantile(fz, .975), slz = mean(slz), slz_lo = quantile(slz, .025), slz_hi = quantile(slz, .975)), by = .(tissue, r)]
for (i in seq_len(nrow(SR))) { x <- SR[i]; sec <- "C. Simulation with equal true slopes in both genotype groups (no cis-by-trans effect)"
  add(sec, sprintf("%s, true co-variation r = %+.2f: correlation-difference statistic", x$tissue, x$r), f4(x$fz), f4(x$fz_lo), f4(x$fz_hi), NA, "50 replicates", "49")
  add(sec, sprintf("%s, true co-variation r = %+.2f: slope-difference statistic", x$tissue, x$r), f4(x$slz), f4(x$slz_lo), f4(x$slz_hi), NA, "50 replicates", "49") }
lg <- readLines("50_h5_bias_expectation.log"); k <- grep("expected_H5 observed_H5_here", lg)
EH <- fread(text = sub("^ *[0-9]+: *", "", lg[(k + 2):(k + 4)]), header = FALSE, col.names = c("sb", "bias_tumour", "bias_normal", "expected_H5", "observed_H5_here"))
for (i in seq_len(nrow(EH))) add("D. Expected tumor-minus-normal difference without cis-by-trans, calibrated to the observed co-variation (tumor +0.004, normal -0.003)",
  sprintf("True SD of |a| = %.2f: expected vs observed (same 30-factor subsample)", EH$sb[i]), sprintf("%s vs %s", f4(EH$expected_H5[i]), f4(EH$observed_H5_here[i])), NA, NA, NA, "", "50")
TS <- fread("prostate_tf_set/set_vs_matched_tf.csv")[stat == "delta"]
for (i in seq_len(nrow(TS))) { x <- TS[i]
  add(sprintf("E. Prostate factor sets versus expression-matched random sets (%s)", c(r0 = "depth only; primary", r_bp = "+ BayesPrism epithelial fraction")[[x$model]]),
      sprintf("%s (%s): tumor minus normal; matched null mean +- SD", c(LIN = "Lineage", ONC = "Driver")[[x$set]], x$members),
      sprintf("%s; null %s +- %.4f", f4(x$obs), f4(x$matched_null_mean), x$matched_null_sd), NA, NA, signif(x$P, 2), "", "52") }
for (tis in c("T", "N")) { z <- fread(sprintf("cistrans_validate_28f/random_nonTF_draws_%s.csv", tis))
  add("F. Genotype x factor interaction-eQTL inflation at GTEx prostate lead variants versus expression-matched non-factor genes",
      sprintf("%s: lambda, transcription factors vs random non-factor genes (mean +- SD; range)", c(T = "Tumor", N = "Normal")[[tis]]),
      sprintf("%.3f vs %.3f +- %.3f (%.3f-%.3f)", z$tf_lambda[1], mean(z$lambda), sd(z$lambda), min(z$lambda), max(z$lambda)), NA, NA,
      signif((sum(z$lambda >= z$tf_lambda[1]) + 1) / (nrow(z) + 1), 2), sprintf("%d draws", nrow(z)), "28f") }
T17 <- rbindlist(rows); fwrite(T17, file.path(O, "TableS17_cis_by_trans_slope.csv"), na = "")
cat(sprintf("S17: %d rows; primary delta %s, permutation P %s\n", nrow(T17), T17[Metric == "Tumor minus normal" & grepl("primary", Section), Value],
            T17[grepl("^Tumor minus normal: whole", Metric) & grepl("primary", Section), P]))
cat("DONE\n")
