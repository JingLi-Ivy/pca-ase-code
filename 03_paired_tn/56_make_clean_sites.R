#!/usr/bin/env Rscript
# 56: clean paired_balanced_sites for 16-bit count overflow (found 2026-10-06; 37b/55):
#  (1) drop rows with any negative count or depth >= 32767 (wrapped / saturated);
#  (2) drop "wrap-risk" genes: any single allele count >= 15000 in any patient (true counts > 65535 wrap to plausible positives).
# Output: delta_ase_clean/paired_balanced_sites.csv.gz (+ list of dropped genes). Downstream *_clean scripts also cap weights at 100.
suppressPackageStartupMessages(library(data.table))
B <- fread("delta_ase/paired_balanced_sites.csv.gz"); n0 <- nrow(B)
ovf <- B[, rnT < 0 | anT < 0 | rnN < 0 | anN < 0 | rnT + anT >= 32767 | rnN + anN >= 32767]
hi <- B[!ovf][pmax(rnT, anT, rnN, anN) >= 15000, unique(eg)]
fwrite(data.table(eg = hi), "delta_ase_clean/wrap_risk_genes.csv")
B <- B[!ovf & !eg %chin% hi]
fwrite(B, "delta_ase_clean/paired_balanced_sites.csv.gz")
cat(sprintf("rows %d -> %d ; overflow rows %d ; wrap-risk genes %d (%s)\n", n0, nrow(B), sum(ovf), length(hi), paste(hi, collapse = ",")))
