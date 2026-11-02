#!/usr/bin/env Rscript
# 68: GO BP recheck for v10 (2026-10-10). Steps 1-4 copied verbatim from figure_freeze/Fig1/code/03_GO_KEGG_enrichment.R (same background and topGO settings);
# only change: BH FDR over ALL tested GO terms (03 applied BH to the 100 terms returned by GenTable(topNodes = 100), underestimating FDR).
# GO/KEGG Enrichment Analysis for PCa-ASE Genes
# Using correct background: all genes with heterozygous SNP coverage in RNA-seq
# FDR threshold: 0.1

.libPaths(c("/path/to/data_root/ASE/code/v8_round2/Rlib_topgo", .libPaths()))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(topGO))
suppressMessages(library(org.Hs.eg.db))

cat("\n==========================================\n")
cat("GO/KEGG Enrichment Analysis for PCa-ASE Genes\n")
cat("==========================================\n\n")

bitr <- function(geneID, fromType, toType, OrgDb) { r <- suppressMessages(AnnotationDbi::select(OrgDb, keys = unique(geneID), keytype = fromType, columns = toType)); unique(r[!is.na(r[[toType]]), c(fromType, toType)]) }

setwd("/path/to/data_root/ASE")

pca_ase_file <- "code/PCa_ASE_genes.csv"
if (!file.exists(pca_ase_file)) {
  stop("Error: PCa-ASE gene list not found: ", pca_ase_file)
}

pca_ase_genes <- read.csv(pca_ase_file, stringsAsFactors = FALSE)
cat(sprintf("PCa-ASE genes: %d\n", nrow(pca_ase_genes)))

pca_ase_gene_symbols <- pca_ase_genes %>%
  filter(!grepl("^ENSG", gene_symbol)) %>%
  pull(gene_symbol) %>%
  unique()

cat(sprintf("PCa-ASE genes with valid symbols: %d\n\n", length(pca_ase_gene_symbols)))

cat("Step 2: Collecting background genes from counts_for_R files...\n")

counts_files <- list.files(pattern = "^counts_for_R_pipeline_.*_with_GQ_geneBody\\.csv$", 
                           full.names = TRUE)
cat(sprintf("Found %d counts files\n", length(counts_files)))

all_background_genes <- c()

for (file in counts_files) {
  tryCatch({
    df <- read.csv(file, stringsAsFactors = FALSE)
    genes <- unique(df$gene_id)
    all_background_genes <- c(all_background_genes, genes)
  }, error = function(e) {
    cat(sprintf("Warning: Failed to read %s\n", basename(file)))
  })
}

all_background_genes <- unique(all_background_genes)
cat(sprintf("Total unique background genes (Ensembl ID): %d\n\n", length(all_background_genes)))

cat("Step 3: Converting gene IDs to Entrez ID...\n")

pca_entrez <- bitr(pca_ase_gene_symbols, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

cat(sprintf("PCa-ASE genes mapped to Entrez ID: %d / %d (%.1f%%)\n", 
            nrow(pca_entrez), 
            length(pca_ase_gene_symbols),
            nrow(pca_entrez) / length(pca_ase_gene_symbols) * 100))

all_background_clean <- gsub("\\..*", "", all_background_genes)
all_background_clean <- unique(all_background_clean)

cat("Converting background genes (Ensembl ID -> Entrez ID)...\n")
universe_entrez_df <- bitr(all_background_clean,
                           fromType = "ENSEMBL",
                           toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)

universe_entrez <- unique(universe_entrez_df$ENTREZID)

cat(sprintf("Background genes mapped to Entrez ID: %d / %d (%.1f%%)\n", 
            length(universe_entrez),
            length(all_background_clean),
            length(universe_entrez) / length(all_background_clean) * 100))
cat(sprintf("Final background gene count (Entrez ID): %d\n\n", length(universe_entrez)))

cat("Step 4: GO Enrichment Analysis (Biological Process)...\n")

geneList <- factor(as.integer(universe_entrez %in% pca_entrez$ENTREZID))
names(geneList) <- universe_entrez

n_interest <- sum(geneList == "1")
n_background <- sum(geneList == "0")
n_total <- length(geneList)

cat(sprintf("  Genes of interest: %d\n", n_interest))
cat(sprintf("  Background genes: %d\n", n_background))
cat(sprintf("  Total genes: %d\n\n", n_total))

GOdata_BP <- new("topGOdata",
                 ontology = "BP",
                 allGenes = geneList,
                 annot = annFUN.org,
                 mapping = "org.Hs.eg.db",
                 ID = "entrez",
                 nodeSize = 5)

result_fisher_BP <- runTest(GOdata_BP, 
                            algorithm = "classic", 
                            statistic = "fisher")

pv <- score(result_fisher_BP); n_tested <- length(pv)
tab <- GenTable(GOdata_BP, Fisher = result_fisher_BP, orderBy = "Fisher", topNodes = n_tested, numChar = 200)
tab$p <- pv[tab$GO.ID]; tab$FDR_all <- p.adjust(tab$p, "BH"); tab$FDR_top100 <- NA
tab$FDR_top100[1:100] <- p.adjust(tab$p[1:100], "BH"); tab$ratio <- tab$Significant / tab$Expected; tab$rank_p <- rank(tab$p, ties.method = "min")
write.csv(tab, "code/v8_round2/supp_data/GO_BP_topGO_classic_allterms_FDR.csv", row.names = FALSE)
cat(sprintf("\nGO BP terms tested (nodeSize 5): %d ; FDR_all < 0.05: %d ; < 0.10: %d (top-100 BH gave %d < 0.10)\n", n_tested, sum(tab$FDR_all < 0.05), sum(tab$FDR_all < 0.10), sum(tab$FDR_top100 < 0.10, na.rm = TRUE)))
for (g in c("GO:0070972", "GO:0048193", "GO:0016192")) { x <- tab[tab$GO.ID == g, ]
  cat(sprintf("%s %s: rank by P %d ; ratio %.2f (%d/%.1f) ; P %.2g ; FDR(all) %.2g ; FDR(top100) %.2g\n", g, x$Term, x$rank_p, x$ratio, x$Significant, x$Expected, x$p, x$FDR_all, x$FDR_top100)) }
cat("\nTerms with FDR(all) < 0.10, ordered by P:\n"); print(tab[tab$FDR_all < 0.10, c("GO.ID", "Term", "Annotated", "Significant", "ratio", "p", "FDR_all")], row.names = FALSE)
cat("DONE\n")
