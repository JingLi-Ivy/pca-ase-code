#!/usr/bin/env Rscript
# 62b: S6 input (ASE cohort, gnomAD AF<0.01): consequence counts once per sample-gene-consequence (fig3A rules) and Venn sets.
# Derived from 62a: can the hard-coded Venn numbers of v9 Fig 3C (SCNA 666, missense 826, overlap 327; Fig3C_venn.R) be reproduced?
# Counts PCa-ASE genes with a somatic missense variant in the per-sample VEP files (fig3A input) under several definitions.
suppressPackageStartupMessages(library(data.table))
VEP <- "/path/to/data_root/wasp/VEP_annotation"
pca <- fread("/path/to/data_root/ASE/code/PCa_ASE_genes.csv")[, g := sub("[.].*", "", gene_id)]
scna <- fread("/path/to/data_root/ASE/code/fig2_results/pca_ase_cnv_correlations.csv")[is.finite(fdr) & fdr < 0.1, gene_id_clean]
coh <- sub("_ase_results[.]csv$", "", list.files("/path/to/data_root/ASE", pattern = "^T[0-9]+_ase_results[.]csv$"))
fs <- list.files(VEP, pattern = "_vep.txt$", full.names = TRUE)
X <- rbindlist(lapply(fs, function(f) { d <- tryCatch(fread(f, skip = "#Uploaded_variation", sep = "\t", na.strings = c("-", ""), quote = "", fill = TRUE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  af <- grep("gnomAD.*_AF$|^AF$", names(d), value = TRUE); mx <- if (length(af)) apply(d[, ..af], 1, function(x) suppressWarnings(max(as.numeric(x), na.rm = TRUE))) else NA
  data.table(sample = paste0("T", sub("_vs_.*", "", basename(f))), var = d[[1]], gene = d$Gene, cons = d$Consequence, impact = if ("IMPACT" %in% names(d)) d$IMPACT else NA, maxaf = mx) }))
X <- X[gene %in% pca$g]
cat("VEP files", length(fs), "; samples", uniqueN(X$sample), "; in ASE cohort", sum(unique(X$sample) %in% coh), "\n")
mis <- X[grepl("missense_variant", cons)]
for (lab in c("all samples, any AF", "all samples, AF<0.01", "ASE cohort, any AF", "ASE cohort, AF<0.01")) {
  m <- mis; if (grepl("AF<", lab)) m <- m[!is.finite(maxaf) | maxaf < 0.01]; if (grepl("cohort", lab)) m <- m[sample %in% coh]
  g <- unique(m$gene); cat(sprintf("%-26s missense PCa-ASE genes %4d ; overlap with 666 SCNA genes %4d ; variants %d\n", lab, length(g), sum(scna %in% g), uniqueN(m[, paste(sample, var)]))) }

pri <- c("stop_gained","frameshift_variant","splice_acceptor_variant","splice_donor_variant","start_lost","stop_lost","inframe_deletion","inframe_insertion","missense_variant")
lab <- c("Stop gained","Frameshift","Splice acceptor","Splice donor","Start lost","Stop lost","Inframe deletion","Inframe insertion","Missense")
C <- X[sample %in% coh & (!is.finite(maxaf) | maxaf < 0.01)]
C[, top := sapply(strsplit(cons, ","), function(v) { h <- match(v, pri); if (all(is.na(h))) NA_character_ else pri[min(h, na.rm = TRUE)] })]
C <- unique(C[!is.na(top), .(sample, gene, top)])
cnt <- C[, .N, by = top][, label := lab[match(top, pri)]][order(-N)]
fwrite(cnt, "supp_data/S6A_consequence_counts_ASEcohort_AF001.csv"); print(cnt)
mg <- unique(C[top == "missense_variant", gene]); fwrite(data.table(set = c("SCNA","missense","overlap"), n = c(length(scna), length(mg), sum(scna %in% mg))), "supp_data/S6C_venn_ASEcohort_AF001.csv")
