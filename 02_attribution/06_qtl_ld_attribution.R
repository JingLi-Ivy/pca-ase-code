#!/usr/bin/env Rscript
# A10: eQTL/sQTL attribution sensitivity: extends the exact-position match between ASE and QTL SNPs to LD (r2>=0.8) and to the gene level
suppressPackageStartupMessages(library(data.table)); setDTthreads(8)
O <- "/path/to/data_root/ASE/code/v8_round2/qtl_ld"; dir.create(O, showWarnings = FALSE); setwd(O)
D <- "/path/to/data_root/CPGEA/QTLs/CPGEA_eQTL_and_sQTL"
strip <- function(x) sub("\\..*", "", sub("\\|.*", "", x))
pca <- fread("/path/to/data_root/ASE/code/PCa_ASE_genes.csv"); G <- unique(strip(pca$gene_id))
cat("PCa-ASE genes:", length(G), "\n")

x <- fread("/path/to/data_root/ASE/code/ASE_SNPs_with_QTL_annotation_exact_match.csv",
           select = c("start","ase_status","gene_id","chr","is_any_eqtl","is_any_sqtl"))
x[, gene := strip(gene_id)]
x <- x[gene %in% G & ase_status == "Significant_ASE"]
exact_orig <- x[, .(exact_orig_e = any(is_any_eqtl), exact_orig_s = any(is_any_sqtl)), by = gene]
cat(sprintf("Original definition reproduced (identical position, any gene): eQTL %d, sQTL %d\n", exact_orig[exact_orig_e == TRUE, .N], exact_orig[exact_orig_s == TRUE, .N]))
ase <- unique(x[, .(gene, chr = as.character(chr), pos = start)])
cat("sig ASE SNP-gene pairs:", nrow(ase), " unique positions:", uniqueN(ase[, .(chr, pos)]), "\n")

rdE <- function(f) { d <- fread(file.path(D, f), select = c("gene","Chr","Pos")); d[, .(gene = strip(gene), chr = sub("chr", "", Chr), pos = Pos)] }
rdS <- function(f) { d <- fread(file.path(D, f), select = c("grp_id","var_chr","var_from")); d[, .(gene = strip(grp_id), chr = as.character(var_chr), pos = var_from)] }
eq <- unique(rbind(rdE("CPGEA_cis_tumor_eQTLs_anno.txt"), rdE("CPGEA_cis_normal_eQTLs_anno.txt")))[gene %in% G]
sq <- unique(rbind(rdS("CPGEA_cis_tumor_sQTLs_anno.txt"), rdS("CPGEA_cis_normal_sQTLs_anno.txt")))[gene %in% G]
cat("QTL variant-gene pairs in PCa-ASE genes: eQTL", nrow(eq), " sQTL", nrow(sq), "\n")

pos_all <- unique(rbind(ase[, .(chr, pos)], eq[, .(chr, pos)], sq[, .(chr, pos)]))
fwrite(pos_all[, .(chr, pos, pos, "x")], "extract.range", sep = "\t", col.names = FALSE)
fwrite(unique(ase[, .(id = paste0(chr, ":", pos))]), "ase_ids.txt", col.names = FALSE)
cmd <- paste("plink2 --bfile /path/to/data_root/CPGEA/CPGEA_genotype/combined.134.qc --extract bed1 extract.range",
             "--set-all-var-ids @:# --rm-dup force-first --make-bed --out sub --memory 16000 --threads 8 > /dev/null",
             "&& plink2 --bfile sub --r2-unphased --ld-snp-list ase_ids.txt --ld-window-kb 1000 --ld-window-r2 0.8",
             "--ld-window 999999 --out ld --memory 16000 --threads 8 > /dev/null")
stopifnot(system(cmd) == 0)
ldf <- list.files(".", pattern = "^ld\\.vcor")[1]; ld <- fread(ldf)
setnames(ld, grep("ID_A", names(ld), value = TRUE), "A"); setnames(ld, grep("ID_B", names(ld), value = TRUE), "B")
ld <- unique(rbind(ld[, .(A, B)], ld[, .(A = B, B = A)], data.table(A = unique(ase[, paste0(chr, ":", pos)]), B = unique(ase[, paste0(chr, ":", pos)]))))
cat("LD pairs (r2>=0.8, including self):", nrow(ld), "\n")

attr_ld <- function(q) {
  a <- ase[, .(gene, A = paste0(chr, ":", pos))]
  m <- merge(a, ld, by = "A", allow.cartesian = TRUE)
  m <- merge(m, q[, .(gene, B = paste0(chr, ":", pos))], by = c("gene", "B"))
  unique(m$gene)
}
exact_same <- function(q) unique(merge(ase, q, by = c("gene","chr","pos"))$gene)
res <- data.table(
  tier = c("exact position, any gene (manuscript)", "exact position, same gene", "LD r2>=0.8, same gene", "gene has any cis-QTL"),
  eQTL = c(exact_orig[exact_orig_e == TRUE, .N], length(exact_same(eq)), length(attr_ld(eq)), uniqueN(eq$gene)),
  sQTL = c(exact_orig[exact_orig_s == TRUE, .N], length(exact_same(sq)), length(attr_ld(sq)), uniqueN(sq$gene)))
res[, eQTL_pct := round(100 * eQTL / length(G), 1)][, sQTL_pct := round(100 * sQTL / length(G), 1)]
print(res)

s5 <- fread("/path/to/data_root/ASE/code/v8_round2/TableS5_source_attribution.csv")
scna <- s5[is_SCNA == TRUE, gene_id_clean]; meth <- s5[is_Methylation == TRUE, gene_id_clean]
cat(sprintf("\nS5: SCNA %d, methylation %d, S5 eQTL %d, S5 sQTL %d\n", length(scna), length(meth), s5[is_eQTL == TRUE, .N], s5[is_sQTL == TRUE, .N]))
un <- function(e, s) sum(!(G %in% c(scna, meth, e, s)))
u <- data.table(tier = res$tier,
  unattributed = c(sum(!(G %in% c(scna, meth, s5[is_eQTL == TRUE, gene_id_clean], s5[is_sQTL == TRUE, gene_id_clean]))),
                   un(exact_same(eq), exact_same(sq)), un(attr_ld(eq), attr_ld(sq)), un(unique(eq$gene), unique(sq$gene))))
u[, pct := round(100 * unattributed / length(G), 1)]
print(u)
fwrite(res, "qtl_attribution_tiers.csv"); fwrite(u, "unattributed_by_tier.csv")
fwrite(data.table(gene = G, eqtl_exact = G %in% exact_same(eq), eqtl_ld = G %in% attr_ld(eq), eqtl_gene = G %in% eq$gene,
                  sqtl_exact = G %in% exact_same(sq), sqtl_ld = G %in% attr_ld(sq), sqtl_gene = G %in% sq$gene), "gene_qtl_tiers.csv")
