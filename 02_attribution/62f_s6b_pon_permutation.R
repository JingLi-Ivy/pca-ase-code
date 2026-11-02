#!/usr/bin/env Rscript
# 62f: S6B on the same mutation set as S6A/S6C (62d: PoN-filtered, ANNOVAR, 118 ASE-cohort tumors), replacing the VEP/no-PoN input of 13/30.
# Categories as in 13 (one cell per tumor x PCa-ASE gene; ASE = any SNP fdr<0.1 & aRNA>0.1; SCNA = FACETS cnv_type != NEUTR).
# Missense (nonsynonymous_SNV) and synonymous (synonymous_SNV, negative control) rates per category; 1000 permutations as Sen 2022:
# within gene (primary, as Sen; controls gene length / mutability) and within tumor (secondary; controls TMB). One-sided P = mean(null >= obs), floor 1/1000.
# Sensitivity: excluding T438/T502 (PoN-filtered call outliers). Gene symbols: first ANNOVAR gene, as 62e.
suppressPackageStartupMessages(library(data.table))
set.seed(123); N_PERM <- 1000
V <- "/path/to/data_root/ASE/code/v8_round2"; P <- "/path/to/data_root/ASE"
pca <- fread(file.path(P, "code/PCa_ASE_genes.csv"))[, g := sub("[.].*", "", gene_id)]
load(file.path(P, "ase_data.RData")); setDT(ase_data)
if (!"aRNA" %in% names(ase_data)) ase_data[, aRNA := abs(a)]
A <- ase_data[sample_type == "Tumor"][, `:=`(sample = sub("^(T[0-9]+).*", "\\1", sample_id), g = sub("[.].*", "", gene_id))][g %in% pca$g]
A <- A[, .(has_ASE = any(fdr < 0.1 & aRNA > 0.1, na.rm = TRUE)), by = .(sample, g)]
fc <- fread("/path/to/data_root/ASE/derived_inputs/gene_cnv_scores_complete.csv")[, g := sub("[.].*", "", gene_id)]
sc <- unique(fc[cnv_type != "NEUTR", .(sample = Sample, g)])[, has_SCNA := TRUE]
D <- merge(A, sc, by = c("sample", "g"), all.x = TRUE)[is.na(has_SCNA), has_SCNA := FALSE]
D[, category := fifelse(!has_ASE, "non_ASE", fifelse(has_SCNA, "ASE_SCNA", "ASE_non_SCNA"))]
X0 <- fread(file.path(V, "supp_data/S6_pon_exonic_variants.tsv"))[, gene := sub("(\\\\x3b|[,;]).*", "", gene)]
cat(sprintf("ASE tumors %d ; PoN tumors %d ; in both %d\n", uniqueN(D$sample), uniqueN(X0$sample), length(intersect(D$sample, X0$sample))))
D <- D[sample %in% X0$sample]
X <- merge(X0, pca[, .(gene = gene_symbol, g)], by = "gene")
for (cl in c("nonsynonymous_SNV", "synonymous_SNV")) { m <- unique(X[exonic_func == cl, .(sample, g)])
  cat(sprintf("%s: %d tumor-gene events in PCa-ASE genes (symbol-matched; %d of %d PCa-ASE symbols ever seen)\n", cl, nrow(m), uniqueN(X$gene), nrow(pca)))
  D[, (cl) := FALSE]; D[m, on = .(sample, g), (cl) := TRUE] }
setnames(D, c("nonsynonymous_SNV", "synonymous_SNV"), c("mis", "syn"))
res <- list(); nulls <- list()
for (excl in list(character(0), c("T438", "T502"))) {
  tag <- if (length(excl)) "excl_T438_T502" else "all"; d <- copy(D[!sample %in% excl])
  obs <- d[, .(n = .N, mis = sum(mis), syn = sum(syn), mis_rate = mean(mis), syn_rate = mean(syn)), by = category]
  print(cbind(set = tag, obs)); o <- obs[category == "ASE_non_SCNA"]
  for (by in c("g", "sample")) {
    nv <- t(vapply(seq_len(N_PERM), function(i) { d[, cp := sample(category), by = by]
      d[cp == "ASE_non_SCNA", c(mean(mis), mean(syn))] }, numeric(2)))
    sch <- if (by == "g") "Within gene" else "Within tumor"
    nulls[[paste(tag, by)]] <- data.table(set = tag, scheme = sch, mis = nv[, 1], syn = nv[, 2])
    res[[paste(tag, by)]] <- data.table(set = tag, scheme = sch, n_cells = o$n, mis_n = o$mis, syn_n = o$syn,
      mis_obs = o$mis_rate, mis_null = mean(nv[, 1]), mis_fold = o$mis_rate / mean(nv[, 1]), mis_P = max(mean(nv[, 1] >= o$mis_rate), 1 / N_PERM),
      syn_obs = o$syn_rate, syn_null = mean(nv[, 2]), syn_fold = o$syn_rate / mean(nv[, 2]), syn_P = max(mean(nv[, 2] >= o$syn_rate), 1 / N_PERM))
  }
}
R <- rbindlist(res); print(R)
fwrite(R, file.path(V, "supp_data/S6B_pon_permutation_summary.csv"))
fwrite(rbindlist(nulls), file.path(V, "supp_data/S6B_pon_permutation_null.csv"))
fwrite(D[, .(sample, g, category, mis, syn)], file.path(V, "supp_data/S6B_pon_per_cell.csv"))
