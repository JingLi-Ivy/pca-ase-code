#!/usr/bin/env Rscript
# 61a: compact plot data for Supplemental Fig S1 (ASE calling QC), recomputed from the SNP-level ASE table with the exact
# procedure and seed of 07_ASE_QC_panel.R (panels A-C; panel D uses 07's saved ase_qc_results/ASE_depth_sensitivity.csv).
#  A depth histogram (counts per depth 1..200, plus >200)   B reference fraction at background sites (depth >= 10), 100 bins
#  C null QQ: beta-binomial (overdispersion by ML on background sites) and binomial two-sided P, 30,000 background sites
suppressPackageStartupMessages(library(data.table))
set.seed(42)
O <- "/path/to/data_root/ASE/code/v8_round2/supp_data"; dir.create(O, showWarnings = FALSE)
snp <- fread("/path/to/data_root/ASE/all_T_samples_SNPs_with_ASE_annotation.csv",
             select = c("ref.matches","alt.matches","N","ref_ratio","alt_ratio","pval","fdr","a","ase_status","gene_id","sample_id"))
snp[, depth := N]; snp[is.na(depth) | depth == 0, depth := ref.matches + alt.matches]
snp[, is_sig := (!is.na(fdr) & fdr < 0.10)]
bg <- snp[is_sig == FALSE & is.finite(ref_ratio) & depth >= 10]
fwrite(snp[is.finite(depth), .N, by = .(depth = pmin(depth, 201))][order(depth)], file.path(O, "S1A_depth_hist.csv"))
fwrite(bg[, .N, by = .(bin = pmin(floor(ref_ratio * 100), 99))][order(bin)], file.path(O, "S1B_refratio_hist.csv"))
fwrite(data.table(n_obs = nrow(snp), n_samples = uniqueN(snp$sample_id), median_depth = median(snp$depth, na.rm = TRUE),
                  bg_n = nrow(bg), bg_median = median(bg$ref_ratio), bg_mean = mean(bg$ref_ratio)), file.path(O, "S1_summary.csv"))
bgc <- bg[is.finite(ref.matches) & depth >= 10]
if (nrow(bgc) > 3e4) bgc <- bgc[sample(.N, 3e4)]
kk <- bgc$ref.matches; nn <- bgc$depth
ei <- sample(length(kk), min(20000, length(kk)))
nll <- function(lr){ rho <- plogis(lr); a <- (1-rho)/(2*rho)
  -sum(lchoose(nn[ei], kk[ei]) + lbeta(kk[ei]+a, nn[ei]-kk[ei]+a) - lbeta(a, a)) }
rho <- plogis(optimize(nll, c(-12, 2))$minimum); a0 <- (1-rho)/(2*rho)
bb_p <- function(n, k, a){ j <- 0:n; l <- lchoose(n,j)+lbeta(j+a,n-j+a)-lbeta(a,a); min(1, sum(exp(l[l <= l[k+1]+1e-9]))) }
pbb  <- mapply(bb_p, nn, kk, MoreArgs = list(a = a0))
pbin <- pmin(1, 2*pmin(pbinom(kk, nn, 0.5), pbinom(kk-1, nn, 0.5, lower.tail = FALSE)))
lam  <- function(p){ p <- pmin(pmax(p,1e-300),1); median(qchisq(1-p,1))/qchisq(0.5,1) }
fwrite(data.table(p_betabinom = pbb, p_binom = pbin), file.path(O, "S1C_null_p.csv"))
fwrite(data.table(rho = rho, lambda_betabinom = lam(pbb), lambda_binom = lam(pbin), n = length(pbb)), file.path(O, "S1C_summary.csv"))
cat(sprintf("rho %.4g ; lambda bb %.3f binom %.3f\n", rho, lam(pbb), lam(pbin)))
