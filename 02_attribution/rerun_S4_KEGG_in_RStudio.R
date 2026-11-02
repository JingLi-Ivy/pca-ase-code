# =============================================================================
# S4: KEGG enrichment of the 620 deletion-driven ASE genes (requires clusterProfiler)
# Purpose: reproduce the autophagy P value with clusterProfiler::enrichKEGG
# Data: pca_ase_cnv_correlations.csv (/path/to/data_root/ASE/code/fig2_results/)
#       when running locally, copy the csv and edit the CORR path
# =============================================================================
suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
})

CORR <- "/path/to/data_root/ASE/code/fig2_results/pca_ase_cnv_correlations.csv"
cc <- read.csv(CORR, stringsAsFactors = FALSE)

del <- subset(cc, is.finite(fdr) & fdr < 0.1 & rho < 0)
cat(sprintf("Deletion-driven genes: %d (expected 620)\n", nrow(del)))

ent <- bitr(unique(del$gene_symbol), fromType = "SYMBOL",
            toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat(sprintf("mapped to Entrez: %d\n", nrow(ent)))

ek <- enrichKEGG(gene = ent$ENTREZID, organism = "hsa",
                 pvalueCutoff = 1, qvalueCutoff = 1)
ek_df <- as.data.frame(ek)
cat("\n=== Top 10 with the default background ===\n")
print(head(ek_df[, c("ID","Description","pvalue","p.adjust","Count")], 10))
cat("\n--- autophagy hsa04140 ---\n")
print(ek_df[ek_df$ID == "hsa04140", c("ID","Description","pvalue","p.adjust","Count")])
write.csv(ek_df, "S4_KEGG_deletion_driven_clusterProfiler_defaultBG.csv", row.names = FALSE)

uni <- bitr(unique(cc$gene_symbol), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
ek2 <- enrichKEGG(gene = ent$ENTREZID, universe = uni, organism = "hsa",
                  pvalueCutoff = 1, qvalueCutoff = 1)
ek2_df <- as.data.frame(ek2)
cat("\n=== Background = 1796 evaluable genes: autophagy ===\n")
print(ek2_df[ek2_df$ID == "hsa04140", c("ID","Description","pvalue","p.adjust","Count")])
write.csv(ek2_df, "S4_KEGG_deletion_driven_clusterProfiler_evalBG.csv", row.names = FALSE)

cat("\nCompare the autophagy P value under both universes.\n")
