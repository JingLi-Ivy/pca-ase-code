#!/usr/bin/env Rscript
# 54: tumour-normal sample identity check. Sites in paired_balanced_sites are germline heterozygous (normal WGS genotype) and
# tumour balanced diploid. If T and N come from the same person, tumour RNA at deep sites should show both alleles about as often
# as normal RNA. A swapped/contaminated pair shows many tumour sites carrying only one allele (expected ~50% for an unrelated person).
# Statistic per patient: fraction of sites (depth >= 30 in both tissues) where a tissue shows <= 1 read of the minor allele.
suppressPackageStartupMessages(library(data.table)); setDTthreads(4)
B <- fread("/path/to/data_root/ASE/code/v8_round2/delta_ase/paired_balanced_sites.csv.gz")
B <- B[rnT >= 0 & anT >= 0 & rnN >= 0 & anN >= 0 & rnT + anT < 32767 & rnN + anN < 32767 & rnT + anT >= 30 & rnN + anN >= 30]
P <- B[, .(sites = .N, monoT = mean(pmin(rnT, anT) <= 1), monoN = mean(pmin(rnN, anN) <= 1)), by = pid][, diff := monoT - monoN]
setorder(P, -diff); fwrite(P, "delta_ase/tn_identity_check.csv")
cat(sprintf("patients %d ; median monoT %.3f monoN %.3f ; diff median %.3f, IQR %.3f\n", nrow(P), median(P$monoT), median(P$monoN), median(P$diff), IQR(P$diff)))
print(head(P, 12)); print(P[pid %in% c(438, 502)])
