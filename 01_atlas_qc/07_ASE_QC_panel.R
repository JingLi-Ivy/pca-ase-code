#!/usr/bin/env Rscript
# ============================================================
# ASE calling QC / reference-bias control panel (Supp Fig S2)
#   Purpose: a QC panel covering the standard ASE concerns
#   (mapping bias / reference bias / null calibration / depth sensitivity /
#    imprinted-gene positive controls).
#
#   Inputs (read-only):
#     - all_T_samples_SNPs_with_ASE_annotation.csv  (SNP-level ASE annotation)
#         columns: ref.matches, alt.matches, N, ref_ratio, alt_ratio,
#             pval, fdr, a, ase_status, gene_id, sample_id, snp_id ...
#     - [optional] imprinted-gene positive controls: MEST_ase_data.csv / H19_ase_data.csv
#   Outputs (written to a new directory):
#     - Fig_ASE_QC_panel.pdf  (six-panel figure A-F)
#     - ASE_QC_summary.csv     (key QC metrics)
#     - ASE_depth_sensitivity.csv
# ============================================================

suppressMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})
set.seed(42)

cat("==========================================\n")
cat("ASE QC / reference-bias control panel\n")
cat("==========================================\n\n")

SNP_ASE_FILE <- "/path/to/data_root/ASE/all_T_samples_SNPs_with_ASE_annotation.csv"
MEST_FILE    <- "/path/to/data_root/ASE/code/MEST_ase_data.csv"
H19_FILE     <- "/path/to/data_root/ASE/code/H19_ase_data.csv"

OUTPUT_DIR <- "/path/to/data_root/ASE/code/ase_qc_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

ASE_FDR        <- 0.10
DEPTH_CUTOFFS  <- c(5, 10, 15, 20, 30)

cat("Step 1: Reading SNP-level ASE annotation (large file, fread with column selection)...\n")
keep <- c("ref.matches","alt.matches","N","ref_ratio","alt_ratio",
          "pval","fdr","a","ase_status","gene_id","sample_id")
snp <- fread(SNP_ASE_FILE, select = keep, showProgress = TRUE)
cat(sprintf("  SNP records: %s\n", format(nrow(snp), big.mark=",")))
cat(sprintf("  Number of samples: %d\n\n", uniqueN(snp$sample_id)))

snp[, depth := N]
snp[is.na(depth) | depth == 0, depth := ref.matches + alt.matches]
snp[, is_sig := (!is.na(fdr) & fdr < ASE_FDR)]
bg <- snp[is_sig == FALSE & is.finite(ref_ratio) & depth >= 10]

cat("Step 2: Panel A: depth distribution...\n")
pA <- ggplot(snp[is.finite(depth) & depth <= 200],
             aes(x = depth)) +
  geom_histogram(bins = 60, fill = "#4393C3", color = NA) +
  geom_vline(xintercept = 10, linetype = 2, color = "red") +
  labs(x = "Allelic read depth (ref+alt)", y = "No. of SNPs",
       title = "RNA allelic read depth") +
  theme_bw(base_size = 11)

cat("Step 3: Panel B: reference bias...\n")
ref_bias_median <- median(bg$ref_ratio, na.rm = TRUE)
ref_bias_mean   <- mean(bg$ref_ratio, na.rm = TRUE)
pB <- ggplot(bg, aes(x = ref_ratio)) +
  geom_histogram(bins = 60, fill = "#B2DF8A", color = NA) +
  geom_vline(xintercept = 0.5, linetype = 2, color = "black") +
  geom_vline(xintercept = ref_bias_median, linetype = 1, color = "red") +
  annotate("text", x = 0.5, y = Inf, vjust = 2,
           label = sprintf("median=%.3f", ref_bias_median), color = "red", size = 3) +
  labs(x = "Reference-allele fraction (non-sig sites)", y = "No. of SNPs",
       title = "Reference-bias check (after WASP)") +
  theme_bw(base_size = 11)

cat("Step 4: Panel C: null QQ (beta-binomial, modelling overdispersion)...\n")
bgc <- bg[is.finite(ref.matches) & depth >= 10]
if (nrow(bgc) > 3e4) bgc <- bgc[sample(.N, 3e4)]
kk <- bgc$ref.matches; nn <- bgc$depth
ei <- sample(length(kk), min(20000, length(kk)))
nll <- function(lr){ rho <- plogis(lr); a <- (1-rho)/(2*rho)
  -sum(lchoose(nn[ei], kk[ei]) + lbeta(kk[ei]+a, nn[ei]-kk[ei]+a) - lbeta(a, a)) }
rho <- plogis(optimize(nll, c(-12, 2))$minimum); a0 <- (1-rho)/(2*rho)
cat(sprintf("  Overdispersion rho = %.4g\n", rho))
bb_p <- function(n, k, a){ j <- 0:n; l <- lchoose(n,j)+lbeta(j+a,n-j+a)-lbeta(a,a)
  min(1, sum(exp(l[l <= l[k+1]+1e-9]))) }
pbb  <- mapply(bb_p, nn, kk, MoreArgs = list(a = a0))
pbin <- pmin(1, 2*pmin(pbinom(kk, nn, 0.5), pbinom(kk-1, nn, 0.5, lower.tail = FALSE)))
lam  <- function(p){ p <- pmin(pmax(p,1e-300),1); median(qchisq(1-p,1))/qchisq(0.5,1) }
lam_bb <- lam(pbb); lam_bn <- lam(pbin)
cat(sprintf("  lambda: beta-binomial=%.3f  binomial=%.3f\n", lam_bb, lam_bn))
mkqq <- function(p, lab){ p <- sort(p[is.finite(p) & p > 0])
  data.frame(model = lab, exp = -log10(ppoints(length(p))), obs = -log10(p)) }
qq <- rbind(mkqq(pbb,  sprintf("beta-binomial (lambda=%.2f)", lam_bb)),
            mkqq(pbin, sprintf("binomial (lambda=%.2f)", lam_bn)))
pC <- ggplot(qq, aes(x = exp, y = obs, color = model)) +
  geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = 2) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_manual(values = c("#1C7293", "#C0392B")) +
  labs(x = expression(Expected~-log[10](P)),
       y = expression(Observed~-log[10](P)),
       title = "Null QQ: beta-binomial vs binomial", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.42, 0.88), legend.text = element_text(size = 7),
        legend.background = element_rect(fill = NA))

cat("Step 5: Panel D: depth sensitivity...\n")
sens <- rbindlist(lapply(DEPTH_CUTOFFS, function(dc) {
  sub <- snp[depth >= dc]
  data.table(min_depth = dc,
             n_snp_sig  = sum(sub$is_sig, na.rm = TRUE),
             n_gene_sig = uniqueN(sub[is_sig == TRUE, gene_id]))
}))
fwrite(sens, file.path(OUTPUT_DIR, "ASE_depth_sensitivity.csv"))
pD <- ggplot(sens, aes(x = factor(min_depth), y = n_gene_sig, group = 1)) +
  geom_line(color = "#2166AC") + geom_point(size = 2, color = "#2166AC") +
  geom_text(aes(label = n_gene_sig), vjust = -0.6, size = 3) +
  labs(x = "Minimum read-depth threshold", y = "No. of significant ASE genes",
       title = "Read-depth sensitivity") +
  theme_bw(base_size = 11)

cat("Step 6: Panel E: imprinted-gene positive control...\n")
read_imprint <- function(f, label) {
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  if ("ref_ratio" %in% colnames(d))  v <- abs(as.numeric(d[["ref_ratio"]]) - 0.5)
  else if ("a" %in% colnames(d))     v <- abs(as.numeric(d[["a"]]))
  else return(NULL)
  data.table(group = label, val = v)
}
bg_a <- data.table(group = "Background",
                   val = abs(snp$ref_ratio[is.finite(snp$ref_ratio) & snp$depth >= 10] - 0.5))
imp <- rbindlist(list(read_imprint(MEST_FILE, "MEST"),
                      read_imprint(H19_FILE,  "H19"),
                      bg_a[sample(.N, min(.N, 1e5))]), fill = TRUE)
if (nrow(imp) > 0) {
  pE <- ggplot(imp, aes(x = group, y = val, fill = group)) +
    geom_violin(scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.size = 0.3) +
    labs(x = NULL, y = "|ref fraction - 0.5| (allelic imbalance)",
         title = "Imprinted-gene positive control") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8))
} else {
  pE <- ggplot() + theme_void() +
    annotate("text", x=0, y=0, label="Imprinting control file missing, skipped")
}

cat("Step 7: Panel F: ref/alt symmetry...\n")
samp <- snp[is.finite(ref.matches) & is.finite(alt.matches) & depth <= 300]
if (nrow(samp) > 1.5e5) samp <- samp[sample(.N, 1.5e5)]
pF <- ggplot(samp, aes(x = ref.matches, y = alt.matches, color = is_sig)) +
  geom_point(size = 0.25, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2) +
  scale_color_manual(values = c(`FALSE`="grey75", `TRUE`="#D6604D"),
                     labels = c("Non-significant","Significant ASE"), name = NULL) +
  labs(x = "Reference-allele reads", y = "Alternate-allele reads",
       title = "Allelic read symmetry") +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.8, 0.2))

qc_summary <- data.frame(
  metric = c("Total SNPs","Number of samples","Median depth","Reference bias (background median)",
             "Reference bias (background mean)","Significant ASE SNPs","Significant ASE genes",
             "Background sites (QQ)"),
  value = c(nrow(snp), uniqueN(snp$sample_id),
            round(median(snp$depth, na.rm=TRUE),1),
            round(ref_bias_median,4), round(ref_bias_mean,4),
            sum(snp$is_sig, na.rm=TRUE),
            uniqueN(snp[is_sig==TRUE, gene_id]),
            nrow(bgc))
)
write.csv(qc_summary, file.path(OUTPUT_DIR, "ASE_QC_summary.csv"), row.names = FALSE)
cat("\nQC summary:\n"); print(qc_summary)

cat("\nStep 8: Writing PDF panel...\n")
pdf(file.path(OUTPUT_DIR, "Fig_ASE_QC_panel.pdf"), width = 12, height = 8)
if (requireNamespace("gridExtra", quietly = TRUE)) {
  gridExtra::grid.arrange(pA, pB, pC, pD, pE, pF, ncol = 3)
} else {
  for (p in list(pA,pB,pC,pD,pE,pF)) print(p)
}
dev.off()

cat("\n==========================================\n")
cat("Done. Output files:\n")
cat(sprintf("  %s/Fig_ASE_QC_panel.pdf\n", OUTPUT_DIR))
cat(sprintf("  %s/ASE_QC_summary.csv\n", OUTPUT_DIR))
cat(sprintf("  %s/ASE_depth_sensitivity.csv\n", OUTPUT_DIR))
cat("==========================================\n")
