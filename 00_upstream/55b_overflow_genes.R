#!/usr/bin/env Rscript
# 55b: map 16-bit overflow sites (55) to genes; overlap with the 2,378 PCa-ASE genes.
suppressPackageStartupMessages(library(data.table))
o <- fread("count_overflow/overflow_sites.tsv", header = FALSE, col.names = c("s","chr","pos","ref","alt","r","a"))
cat("rows", nrow(o), " negative", o[, sum(r < 0 | a < 0)], " samples", uniqueN(o$s), " sites", uniqueN(o[, paste(chr, pos)]), "\n")
an <- unique(fread("/path/to/data_root/ASE/all_T_samples_SNPs_with_ASE_annotation.csv", select = c("seqtnames","start","gene_id")))
setnames(an, c("chr","pos","gene_id")); an[, chr := as.character(chr)]; o[, chr := as.character(chr)]
g <- merge(o, an, by = c("chr","pos"))[, gene_id := sub("[.].*", "", gene_id)][, .(n_rows = .N, n_neg = sum(r < 0 | a < 0), samples = uniqueN(s)), by = gene_id]
p <- fread("/path/to/data_root/ASE/code/PCa_ASE_genes.csv")[, gid := sub("[.].*", "", gene_id)]
g[, pca_ase := gene_id %in% p$gid]; fwrite(g[order(-n_rows)], "count_overflow/overflow_genes.csv")
cat("overflow genes", nrow(g), "; in 2,378 PCa-ASE:", sum(g$pca_ase), "\n"); print(p[gid %in% g$gene_id, .(rank, gene_symbol, n_sig_Normal, n_sig_Tumor, pct_sig_Tumor)])
