#!/usr/bin/env Rscript
# 00 (2026-10-27): the 2,378-gene PCa-ASE catalogue (Supp. Table S2, code/PCa_ASE_genes.csv), reconstructed because the
#     script that wrote the table on 2025-12-19 was not kept. Reproduces the submitted table exactly (checked below).
# Input: ase_data.RData, saved 2025-12-19 10:45 by code/"View ase_by_expression_quartiles.Rmd" (step 1): all 236
#     {T,N}{ID}_ase_results.csv rows with non-missing pval, fdr and a; ASE = |a|; gene_symbol from a biomaRt query of the
#     Ensembl release current at that date (gene_id where no name was returned). The frozen object is used because that
#     release cannot be queried again.
# Definition (Methods; Sen et al. 2022 neuroblastoma criteria), applied per gene symbol and tissue:
#     n_samples = distinct samples with a test; n_sig = tests with FDR < 0.1; mean_ASE = mean |a| over tests.
#     PCa-ASE gene: n_samples >= 10 in tumor and normal, n_sig >= 3 in tumor and <= 1 in normal.
#     Ten catalogue symbols cover more than one Ensembl ID (e.g. HERC3, U8); their tests are pooled, so a sample can
#     contribute two tests. gene_id = first ID of the symbol in ase_data. Rows ordered by |diff_mean_ASE|;
#     enrichment = (n_sig_T / n_T) / (max(n_sig_N, 1) / n_N).
# Output: OUT/PCa_ASE_genes_reconstructed.csv and OUT/comparison.txt.
suppressPackageStartupMessages(library(data.table))
RDATA <- "/path/to/data_root/ASE/ase_data.RData"
OLD   <- "/path/to/data_root/ASE/code/PCa_ASE_genes.csv"
OUT   <- "/path/to/data_root/ASE/code/v8_round2/catalogue_rederive_20261027"; dir.create(OUT, showWarnings = FALSE)

load(RDATA); d <- as.data.table(ase_data)
g <- d[, .(n = uniqueN(sample_id), n_sig = sum(fdr < 0.1), mean_a = mean(ASE)), by = .(gene_symbol, sample_type)]
w <- dcast(g, gene_symbol ~ sample_type, value.var = c("n", "n_sig", "mean_a"))
w <- w[n_Normal >= 10 & n_Tumor >= 10 & n_sig_Tumor >= 3 & n_sig_Normal <= 1]
w <- merge(w, d[, .(gene_id = gene_id[1]), by = gene_symbol], by = "gene_symbol")
cat_ <- w[, .(gene_id, gene_symbol, n_samples_Normal = n_Normal, n_samples_Tumor = n_Tumor,
              n_sig_Normal, n_sig_Tumor, pct_sig_Normal = 100 * n_sig_Normal / n_Normal, pct_sig_Tumor = 100 * n_sig_Tumor / n_Tumor,
              mean_ASE_Normal = mean_a_Normal, mean_ASE_Tumor = mean_a_Tumor, diff_mean_ASE = mean_a_Tumor - mean_a_Normal,
              enrichment = (n_sig_Tumor / n_Tumor) / (pmax(n_sig_Normal, 1) / n_Normal))]
cat_ <- cat_[order(-abs(diff_mean_ASE))]
cat_ <- cbind(rank = seq_len(nrow(cat_)), cat_)
fwrite(cat_, file.path(OUT, "PCa_ASE_genes_reconstructed.csv"))

old <- fread(OLD)
m <- merge(old, cat_, by = "gene_id", suffixes = c(".old", ".new"))
sink(file.path(OUT, "comparison.txt"), split = TRUE)
cat(sprintf("submitted %d; reconstructed %d; shared gene_id %d\n", nrow(old), nrow(cat_), nrow(m)))
for (v in setdiff(names(old), "gene_id")) {
  a <- m[[paste0(v, ".old")]]; b <- m[[paste0(v, ".new")]]
  cat(sprintf("  %-17s %s\n", v, if (is.numeric(a)) sprintf("max |difference| %.3g", max(abs(a - b))) else sprintf("identical %d / %d", sum(a == b), length(a))))
}
sink()
