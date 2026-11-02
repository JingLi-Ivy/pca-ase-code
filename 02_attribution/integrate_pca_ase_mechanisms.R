#!/usr/bin/env Rscript
# ============================================
# PCa-ASE gene mechanism integration
# Integrates SCNA, eQTL, sQTL, GWAS and methylation sources
# ============================================

library(data.table)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(UpSetR)

cat("============================================\n")
cat("PCa-ASE Genes Mechanism Integration\n")
cat("============================================\n\n")

PCA_ASE_GENES_FILE <- "/path/to/data_root/ASE/code/PCa_ASE_genes.csv"

SCNA_CORR_FILE <- "/path/to/data_root/ASE/code/fig2_results/pca_ase_cnv_correlations.csv"

QTL_ANNOTATION_FILE <- "/path/to/data_root/ASE/code/ASE_SNPs_with_QTL_annotation_exact_match.csv"

GWAS_PCA_FILE <- "/path/to/data_root/ASE/code/PCA_ASE_SNPs_with_GWAS_fixed.csv"

METH_FULL_FILE <- "/path/to/data_root/variants/batch_43/joint_genotyping/pca_ase_methylation_analysis/pca_ase_methylation_analysis_full.txt"
METH_SIG_FILE <- "/path/to/data_root/variants/batch_43/joint_genotyping/pca_ase_methylation_analysis/pca_ase_genes_due_to_methylation.txt"

OUTPUT_DIR <- "/path/to/data_root/ASE/code/mechanism_integration"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("Step 1: Loading PCa-ASE genes...\n")

pca_ase_genes <- fread(PCA_ASE_GENES_FILE)
cat(sprintf("  Total PCa-ASE genes: %d\n", nrow(pca_ase_genes)))

if (!"gene_id_clean" %in% names(pca_ase_genes)) {
  pca_ase_genes[, gene_id_clean := sub("\\..*", "", gene_id)]
}

integration_results <- pca_ase_genes[, .(gene_id, gene_id_clean, gene_symbol)]
integration_results <- unique(integration_results, by = "gene_id_clean")

cat(sprintf("  Unique genes for integration: %d\n\n", nrow(integration_results)))

cat("Step 2: Integrating SCNA-related genes...\n")

if (file.exists(SCNA_CORR_FILE)) {
  scna_corr <- fread(SCNA_CORR_FILE)

  cat("  SCNA file columns:", paste(colnames(scna_corr), collapse = ", "), "\n")

  if("gene_id_clean" %in% colnames(scna_corr) && "fdr" %in% colnames(scna_corr)) {
    scna_genes <- scna_corr[fdr < 0.1]

    if(nrow(scna_genes) > 0) {
      scna_summary <- scna_genes[, .(
        scna_rho = mean(rho, na.rm = TRUE),
        scna_pval = min(pval, na.rm = TRUE),
        scna_fdr = min(fdr, na.rm = TRUE)
      ), by = gene_id_clean]

      scna_summary[, is_SCNA := TRUE]

      cat(sprintf("  SCNA-related genes (FDR<0.1): %d\n", nrow(scna_summary)))

      integration_results <- merge(integration_results, scna_summary, 
                                   by = "gene_id_clean", all.x = TRUE)
      integration_results[is.na(is_SCNA), is_SCNA := FALSE]
    } else {
      cat("  No SCNA-related genes with FDR<0.1\n")
      integration_results[, is_SCNA := FALSE]
    }
  } else {
    cat("  WARNING: Cannot find required columns (gene_id_clean or fdr) in SCNA file\n")
    integration_results[, is_SCNA := FALSE]
  }
} else {
  cat("  WARNING: SCNA file not found\n")
  integration_results[, is_SCNA := FALSE]
}

cat("\nStep 3: Integrating eQTL-related genes...\n")

if (file.exists(QTL_ANNOTATION_FILE)) {
  qtl_data <- fread(QTL_ANNOTATION_FILE, nrows = 1000)
  cat("  QTL file columns:", paste(head(colnames(qtl_data), 10), collapse = ", "), "...\n")

  qtl_data <- fread(QTL_ANNOTATION_FILE)

  ase_col <- if("ase_status" %in% colnames(qtl_data)) "ase_status" else NULL
  eqtl_tumor_col <- if("is_tumor_eqtl" %in% colnames(qtl_data)) "is_tumor_eqtl" else if("tumor_eqtl" %in% colnames(qtl_data)) "tumor_eqtl" else NULL
  eqtl_normal_col <- if("is_normal_eqtl" %in% colnames(qtl_data)) "is_normal_eqtl" else if("normal_eqtl" %in% colnames(qtl_data)) "normal_eqtl" else NULL

  if(!is.null(ase_col) && (!is.null(eqtl_tumor_col) || !is.null(eqtl_normal_col))) {
    if(!is.null(eqtl_tumor_col) && !is.null(eqtl_normal_col)) {
      eqtl_genes <- qtl_data[
        (get(ase_col) == "Significant_ASE" | get(ase_col) == TRUE) & 
        (get(eqtl_tumor_col) == TRUE | get(eqtl_normal_col) == TRUE)
      ]
    } else if(!is.null(eqtl_tumor_col)) {
      eqtl_genes <- qtl_data[
        (get(ase_col) == "Significant_ASE" | get(ase_col) == TRUE) & 
        get(eqtl_tumor_col) == TRUE
      ]
    } else {
      eqtl_genes <- qtl_data[
        (get(ase_col) == "Significant_ASE" | get(ase_col) == TRUE) & 
        get(eqtl_normal_col) == TRUE
      ]
    }

    if(nrow(eqtl_genes) > 0) {
      gene_id_col <- if("gene_id" %in% colnames(eqtl_genes)) "gene_id" else if("gene_name" %in% colnames(eqtl_genes)) "gene_name" else NULL

      if(!is.null(gene_id_col)) {
        eqtl_summary <- eqtl_genes[, .(
          n_eqtl_snps = .N,
          has_tumor_eqtl = if(!is.null(eqtl_tumor_col)) any(get(eqtl_tumor_col) == TRUE, na.rm = TRUE) else FALSE,
          has_normal_eqtl = if(!is.null(eqtl_normal_col)) any(get(eqtl_normal_col) == TRUE, na.rm = TRUE) else FALSE
        ), by = gene_id]

        eqtl_summary[, gene_id_clean := sub("\\..*", "", gene_id)]
        eqtl_summary[, is_eQTL := TRUE]

        cat(sprintf("  eQTL-related genes: %d\n", nrow(eqtl_summary)))
        if(!is.null(eqtl_tumor_col) && !is.null(eqtl_normal_col)) {
          cat(sprintf("    Tumor eQTL only: %d\n", sum(eqtl_summary$has_tumor_eqtl & !eqtl_summary$has_normal_eqtl)))
          cat(sprintf("    Normal eQTL only: %d\n", sum(!eqtl_summary$has_tumor_eqtl & eqtl_summary$has_normal_eqtl)))
          cat(sprintf("    Both: %d\n", sum(eqtl_summary$has_tumor_eqtl & eqtl_summary$has_normal_eqtl)))
        }

        integration_results <- merge(
          integration_results, 
          eqtl_summary[, .(gene_id_clean, n_eqtl_snps, has_tumor_eqtl, has_normal_eqtl, is_eQTL)],
          by = "gene_id_clean", all.x = TRUE
        )
        integration_results[is.na(is_eQTL), is_eQTL := FALSE]
      } else {
        cat("  WARNING: Cannot find gene_id column in QTL file\n")
        integration_results[, is_eQTL := FALSE]
      }
    } else {
      cat("  No eQTL genes found\n")
      integration_results[, is_eQTL := FALSE]
    }
  } else {
    cat("  WARNING: Cannot find required columns in QTL file\n")
    integration_results[, is_eQTL := FALSE]
  }
} else {
  cat("  WARNING: QTL annotation file not found\n")
  integration_results[, is_eQTL := FALSE]
}

cat("\nStep 4: Integrating sQTL-related genes...\n")

if (file.exists(QTL_ANNOTATION_FILE)) {
  qtl_data <- fread(QTL_ANNOTATION_FILE)

  ase_col <- if("ase_status" %in% colnames(qtl_data)) "ase_status" else NULL
  sqtl_tumor_col <- if("is_tumor_sqtl" %in% colnames(qtl_data)) "is_tumor_sqtl" else if("tumor_sqtl" %in% colnames(qtl_data)) "tumor_sqtl" else NULL
  sqtl_normal_col <- if("is_normal_sqtl" %in% colnames(qtl_data)) "is_normal_sqtl" else if("normal_sqtl" %in% colnames(qtl_data)) "normal_sqtl" else NULL

  if(!is.null(ase_col) && (!is.null(sqtl_tumor_col) || !is.null(sqtl_normal_col))) {
    if(!is.null(sqtl_tumor_col) && !is.null(sqtl_normal_col)) {
      sqtl_genes <- qtl_data[
        (get(ase_col) == "Significant_ASE" | get(ase_col) == TRUE) & 
        (get(sqtl_tumor_col) == TRUE | get(sqtl_normal_col) == TRUE)
      ]
    } else if(!is.null(sqtl_tumor_col)) {
      sqtl_genes <- qtl_data[
        (get(ase_col) == "Significant_ASE" | get(ase_col) == TRUE) & 
        get(sqtl_tumor_col) == TRUE
      ]
    } else {
      sqtl_genes <- qtl_data[
        (get(ase_col) == "Significant_ASE" | get(ase_col) == TRUE) & 
        get(sqtl_normal_col) == TRUE
      ]
    }

    if(nrow(sqtl_genes) > 0) {
      gene_id_col <- if("gene_id" %in% colnames(sqtl_genes)) "gene_id" else if("gene_name" %in% colnames(sqtl_genes)) "gene_name" else NULL

      if(!is.null(gene_id_col)) {
        sqtl_summary <- sqtl_genes[, .(
          n_sqtl_snps = .N,
          has_tumor_sqtl = if(!is.null(sqtl_tumor_col)) any(get(sqtl_tumor_col) == TRUE, na.rm = TRUE) else FALSE,
          has_normal_sqtl = if(!is.null(sqtl_normal_col)) any(get(sqtl_normal_col) == TRUE, na.rm = TRUE) else FALSE
        ), by = gene_id]

        sqtl_summary[, gene_id_clean := sub("\\..*", "", gene_id)]
        sqtl_summary[, is_sQTL := TRUE]

        cat(sprintf("  sQTL-related genes: %d\n", nrow(sqtl_summary)))

        integration_results <- merge(
          integration_results,
          sqtl_summary[, .(gene_id_clean, n_sqtl_snps, has_tumor_sqtl, has_normal_sqtl, is_sQTL)],
          by = "gene_id_clean", all.x = TRUE
        )
        integration_results[is.na(is_sQTL), is_sQTL := FALSE]
      } else {
        integration_results[, is_sQTL := FALSE]
      }
    } else {
      cat("  No sQTL genes found\n")
      integration_results[, is_sQTL := FALSE]
    }
  } else {
    integration_results[, is_sQTL := FALSE]
  }
} else {
  integration_results[, is_sQTL := FALSE]
}

cat("\nStep 5: Integrating GWAS-related genes...\n")

if (file.exists(GWAS_PCA_FILE)) {
  gwas_data <- fread(GWAS_PCA_FILE, nrows = 1000)
  cat("  GWAS file columns:", paste(head(colnames(gwas_data), 10), collapse = ", "), "...\n")

  gwas_data <- fread(GWAS_PCA_FILE)

  pval_col <- if("pvalue" %in% colnames(gwas_data)) "pvalue" else if("pval" %in% colnames(gwas_data)) "pval" else if("p.value" %in% colnames(gwas_data)) "p.value" else NULL

  if(!is.null(pval_col)) {
    gwas_sig <- gwas_data[get(pval_col) < 5e-8]

    if (nrow(gwas_sig) > 0) {
      gene_id_col <- if("gene_id" %in% colnames(gwas_sig)) "gene_id" else if("gene_name" %in% colnames(gwas_sig)) "gene_name" else NULL

      if (!is.null(gene_id_col)) {
        gwas_genes <- gwas_sig[, .(
          n_gwas_snps = .N,
          min_gwas_pval = min(get(pval_col), na.rm = TRUE)
        ), by = get(gene_id_col)]

        setnames(gwas_genes, "V1", gene_id_col)
        gwas_genes[, gene_id_clean := sub("\\..*", "", get(gene_id_col))]
        gwas_genes[, is_GWAS := TRUE]

        cat(sprintf("  GWAS-significant genes (p<5e-8): %d\n", nrow(gwas_genes)))

        integration_results <- merge(
          integration_results,
          gwas_genes[, .(gene_id_clean, n_gwas_snps, min_gwas_pval, is_GWAS)],
          by = "gene_id_clean", all.x = TRUE
        )
        integration_results[is.na(is_GWAS), is_GWAS := FALSE]
      } else {
        cat("  Note: gene_id not in GWAS file, matching via QTL file...\n")

        if (file.exists(QTL_ANNOTATION_FILE)) {
          qtl_snp_map <- fread(QTL_ANNOTATION_FILE, select = c("pos_id", "gene_id"))
          qtl_snp_map <- unique(qtl_snp_map[!is.na(gene_id)])
          qtl_snp_map[, gene_id_clean := sub("\\..*", "", gene_id)]

          gwas_sig[, pos_id := paste0(seqtnames, ":", start)]
          gwas_genes <- merge(
            gwas_sig[, .(pos_id, pvalue = get(pval_col))],
            qtl_snp_map[, .(pos_id, gene_id_clean)],
            by = "pos_id"
          )

          if (nrow(gwas_genes) > 0) {
            gwas_summary <- gwas_genes[, .(
              n_gwas_snps = .N,
              min_gwas_pval = min(pvalue, na.rm = TRUE)
            ), by = gene_id_clean]

            gwas_summary[, is_GWAS := TRUE]

            cat(sprintf("  GWAS-significant genes (matched via QTL): %d\n", nrow(gwas_summary)))

            integration_results <- merge(
              integration_results,
              gwas_summary,
              by = "gene_id_clean", all.x = TRUE
            )
            integration_results[is.na(is_GWAS), is_GWAS := FALSE]
          } else {
            cat("  No GWAS genes matched via QTL file\n")
            integration_results[, is_GWAS := FALSE]
          }
        } else {
          cat("  Cannot match GWAS: QTL file not available\n")
          integration_results[, is_GWAS := FALSE]
        }
      }
    } else {
      cat("  No GWAS-significant SNPs found\n")
      integration_results[, is_GWAS := FALSE]
    }
  } else {
    cat("  WARNING: Cannot find pvalue column in GWAS file\n")
    integration_results[, is_GWAS := FALSE]
  }
} else {
  cat("  WARNING: GWAS file not found\n")
  integration_results[, is_GWAS := FALSE]
}

cat("\nStep 6: Integrating Methylation-related genes...\n")

if (file.exists(METH_FULL_FILE)) {
  meth_full <- fread(METH_FULL_FILE)

  cat("  Methylation file columns:", paste(colnames(meth_full), collapse = ", "), "\n")

  fdr_col <- if("fdr" %in% colnames(meth_full)) "fdr" else if("FDR" %in% colnames(meth_full)) "FDR" else NULL

  if(!is.null(fdr_col)) {
    meth_sig <- meth_full[!is.na(get(fdr_col)) & get(fdr_col) < 0.1]

    if (nrow(meth_sig) > 0) {
      meth_genes <- meth_sig[, .(
        gene_id_clean = sub("\\..*", "", gene_id),
        meth_rho = if("rho" %in% colnames(meth_sig)) rho else NA,
        meth_fdr = get(fdr_col),
        n_dmr = if("n_dmr" %in% colnames(meth_sig)) n_dmr else NA,
        dmr_types = if("dmr_types" %in% colnames(meth_sig)) dmr_types else NA
      )]

      meth_genes <- meth_genes[, .(
        meth_rho = mean(meth_rho, na.rm = TRUE),
        meth_fdr = min(meth_fdr, na.rm = TRUE),
        n_dmr = max(n_dmr, na.rm = TRUE),
        dmr_types = paste(unique(dmr_types[!is.na(dmr_types)]), collapse = ";")
      ), by = gene_id_clean]

      meth_genes[, is_Methylation := TRUE]

      cat(sprintf("  Methylation-related genes (FDR<0.1): %d\n", nrow(meth_genes)))

      integration_results <- merge(
        integration_results,
        meth_genes,
        by = "gene_id_clean", all.x = TRUE
      )
      integration_results[is.na(is_Methylation), is_Methylation := FALSE]
    } else {
      cat("  No significant methylation-related genes (FDR<0.1)\n")
      integration_results[, is_Methylation := FALSE]
    }
  } else {
    cat("  WARNING: Cannot find fdr column in methylation file\n")
    integration_results[, is_Methylation := FALSE]
  }
} else {
  cat("  WARNING: Methylation file not found\n")
  integration_results[, is_Methylation := FALSE]
}

cat("\nStep 7: Summarizing mechanisms...\n")

integration_results[, n_mechanisms := 
  as.integer(is_SCNA) + 
  as.integer(is_eQTL) + 
  as.integer(is_sQTL) + 
  as.integer(is_GWAS) + 
  as.integer(is_Methylation)
]

integration_results[, mechanism_combo := paste0(
  ifelse(is_SCNA, "SCNA,", ""),
  ifelse(is_eQTL, "eQTL,", ""),
  ifelse(is_sQTL, "sQTL,", ""),
  ifelse(is_GWAS, "GWAS,", ""),
  ifelse(is_Methylation, "Meth,", "")
)]

integration_results[, mechanism_combo := sub(",$", "", mechanism_combo)]
integration_results[mechanism_combo == "", mechanism_combo := "Unknown"]

integration_results[, primary_mechanism := fifelse(
  is_SCNA, "SCNA",
  fifelse(is_eQTL, "eQTL",
  fifelse(is_sQTL, "sQTL",
  fifelse(is_GWAS, "GWAS",
  fifelse(is_Methylation, "Methylation", "Unknown")))))
]

cat("\n============================================\n")
cat("Summary Statistics\n")
cat("============================================\n\n")

cat("Total PCa-ASE genes:", nrow(integration_results), "\n\n")

cat("Genes by mechanism:\n")
cat(sprintf("  SCNA-related:        %d (%.1f%%)\n", 
            sum(integration_results$is_SCNA), 
            100 * sum(integration_results$is_SCNA) / nrow(integration_results)))
cat(sprintf("  eQTL-related:        %d (%.1f%%)\n",
            sum(integration_results$is_eQTL),
            100 * sum(integration_results$is_eQTL) / nrow(integration_results)))
cat(sprintf("  sQTL-related:        %d (%.1f%%)\n",
            sum(integration_results$is_sQTL),
            100 * sum(integration_results$is_sQTL) / nrow(integration_results)))
cat(sprintf("  GWAS-related:        %d (%.1f%%)\n",
            sum(integration_results$is_GWAS),
            100 * sum(integration_results$is_GWAS) / nrow(integration_results)))
cat(sprintf("  Methylation-related: %d (%.1f%%)\n",
            sum(integration_results$is_Methylation),
            100 * sum(integration_results$is_Methylation) / nrow(integration_results)))
cat(sprintf("  Unknown mechanism:   %d (%.1f%%)\n",
            sum(integration_results$n_mechanisms == 0),
            100 * sum(integration_results$n_mechanisms == 0) / nrow(integration_results)))

cat("\nGenes by number of mechanisms:\n")
print(table(integration_results$n_mechanisms))

cat("\nTop mechanism combinations:\n")
combo_table <- integration_results[, .N, by = mechanism_combo][order(-N)]
print(head(combo_table, 15))

cat("\nStep 8: Generating visualizations...\n")

mechanism_counts <- data.frame(
  Mechanism = c("SCNA", "eQTL", "sQTL", "GWAS", "Methylation", "Unknown"),
  Count = c(
    sum(integration_results$is_SCNA),
    sum(integration_results$is_eQTL),
    sum(integration_results$is_sQTL),
    sum(integration_results$is_GWAS),
    sum(integration_results$is_Methylation),
    sum(integration_results$n_mechanisms == 0)
  )
)

mechanism_counts$Percentage <- 100 * mechanism_counts$Count / nrow(integration_results)
mechanism_counts$Mechanism <- factor(mechanism_counts$Mechanism, 
                                     levels = mechanism_counts$Mechanism[order(-mechanism_counts$Count)])

p1 <- ggplot(mechanism_counts, aes(x = Mechanism, y = Count, fill = Mechanism)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)), 
            vjust = -0.3, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "PCa-ASE Genes by ASE Mechanism",
    subtitle = sprintf("Total genes: %d", nrow(integration_results)),
    x = "Mechanism",
    y = "Number of Genes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ylim(0, max(mechanism_counts$Count) * 1.15)

ggsave(file.path(OUTPUT_DIR, "mechanism_barplot.pdf"), p1, width = 10, height = 6)

upset_data <- integration_results[, .(
  SCNA = as.integer(is_SCNA),
  eQTL = as.integer(is_eQTL),
  sQTL = as.integer(is_sQTL),
  GWAS = as.integer(is_GWAS),
  Methylation = as.integer(is_Methylation)
)]

pdf(file.path(OUTPUT_DIR, "mechanism_upset.pdf"), width = 12, height = 8)
upset(
  as.data.frame(upset_data),
  sets = c("SCNA", "eQTL", "sQTL", "GWAS", "Methylation"),
  order.by = "freq",
  decreasing = TRUE,
  mb.ratio = c(0.6, 0.4),
  text.scale = 1.5,
  main.bar.color = "#4292C6",
  sets.bar.color = "#2171B5",
  matrix.color = "#2171B5"
)
dev.off()

primary_counts <- integration_results[, .N, by = primary_mechanism]
primary_counts[, pct := 100 * N / sum(N)]

p3 <- ggplot(primary_counts, aes(x = "", y = N, fill = primary_mechanism)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  geom_text(aes(label = sprintf("%s\n%d (%.1f%%)", primary_mechanism, N, pct)),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Primary Mechanism Distribution") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(OUTPUT_DIR, "primary_mechanism_pie.pdf"), p3, width = 8, height = 8)

cat("  Visualizations saved\n")

cat("\nStep 9: Saving results...\n")

setorder(integration_results, -n_mechanisms, primary_mechanism, gene_symbol)

fwrite(integration_results, file.path(OUTPUT_DIR, "pca_ase_genes_mechanism_integration.csv"))
cat("  Saved: pca_ase_genes_mechanism_integration.csv\n")

for (mech in c("SCNA", "eQTL", "sQTL", "GWAS", "Methylation")) {
  mech_col <- paste0("is_", mech)
  mech_genes <- integration_results[get(mech_col) == TRUE]
  if (nrow(mech_genes) > 0) {
    fwrite(mech_genes, file.path(OUTPUT_DIR, paste0("genes_", mech, ".csv")))
    cat(sprintf("  Saved: genes_%s.csv (%d genes)\n", mech, nrow(mech_genes)))
  }
}

unknown_genes <- integration_results[n_mechanisms == 0]
if (nrow(unknown_genes) > 0) {
  fwrite(unknown_genes, file.path(OUTPUT_DIR, "genes_unknown_mechanism.csv"))
  cat(sprintf("  Saved: genes_unknown_mechanism.csv (%d genes)\n", nrow(unknown_genes)))
}

multi_mech_genes <- integration_results[n_mechanisms >= 2]
if (nrow(multi_mech_genes) > 0) {
  fwrite(multi_mech_genes, file.path(OUTPUT_DIR, "genes_multiple_mechanisms.csv"))
  cat(sprintf("  Saved: genes_multiple_mechanisms.csv (%d genes)\n", nrow(multi_mech_genes)))
}

summary_stats <- data.table(
  Category = c("Total PCa-ASE genes",
               "SCNA-related", "eQTL-related", "sQTL-related", 
               "GWAS-related", "Methylation-related", "Unknown mechanism",
               "Single mechanism", "Multiple mechanisms (>=2)"),
  Count = c(
    nrow(integration_results),
    sum(integration_results$is_SCNA),
    sum(integration_results$is_eQTL),
    sum(integration_results$is_sQTL),
    sum(integration_results$is_GWAS),
    sum(integration_results$is_Methylation),
    sum(integration_results$n_mechanisms == 0),
    sum(integration_results$n_mechanisms == 1),
    sum(integration_results$n_mechanisms >= 2)
  )
)
summary_stats[, Percentage := sprintf("%.1f%%", 100 * Count / nrow(integration_results))]

fwrite(summary_stats, file.path(OUTPUT_DIR, "mechanism_summary_stats.csv"))
cat("  Saved: mechanism_summary_stats.csv\n")

cat("\nStep 10: Cross-tabulation analysis...\n")

mechanisms <- c("is_SCNA", "is_eQTL", "is_sQTL", "is_GWAS", "is_Methylation")
mech_names <- c("SCNA", "eQTL", "sQTL", "GWAS", "Methylation")

cross_tab <- matrix(0, nrow = 5, ncol = 5)
rownames(cross_tab) <- mech_names
colnames(cross_tab) <- mech_names

for (i in 1:5) {
  for (j in 1:5) {
    cross_tab[i, j] <- sum(integration_results[[mechanisms[i]]] & 
                           integration_results[[mechanisms[j]]])
  }
}

cat("\nMechanism co-occurrence matrix:\n")
print(cross_tab)

cross_tab_dt <- as.data.table(cross_tab, keep.rownames = "Mechanism")
fwrite(cross_tab_dt, file.path(OUTPUT_DIR, "mechanism_cooccurrence_matrix.csv"))

cat("\nFisher's exact test for mechanism associations:\n")

fisher_results <- list()

for (i in 1:(length(mechanisms)-1)) {
  for (j in (i+1):length(mechanisms)) {
    contingency <- table(
      integration_results[[mechanisms[i]]],
      integration_results[[mechanisms[j]]]
    )

    if (all(dim(contingency) == c(2, 2))) {
      fisher_test <- fisher.test(contingency)

      fisher_results[[paste(mech_names[i], mech_names[j], sep = "_vs_")]] <- data.table(
        Mechanism1 = mech_names[i],
        Mechanism2 = mech_names[j],
        Overlap = cross_tab[i, j],
        OR = fisher_test$estimate,
        pvalue = fisher_test$p.value,
        CI_lower = fisher_test$conf.int[1],
        CI_upper = fisher_test$conf.int[2]
      )
    }
  }
}

fisher_dt <- rbindlist(fisher_results)
fisher_dt[, FDR := p.adjust(pvalue, method = "BH")]
fisher_dt <- fisher_dt[order(pvalue)]

cat("\nSignificant associations (FDR < 0.05):\n")
print(fisher_dt[FDR < 0.05])

fwrite(fisher_dt, file.path(OUTPUT_DIR, "mechanism_association_fisher_test.csv"))

cat("\nStep 11: Generating co-occurrence heatmap...\n")

library(pheatmap)

jaccard_matrix <- matrix(0, nrow = 5, ncol = 5)
rownames(jaccard_matrix) <- mech_names
colnames(jaccard_matrix) <- mech_names

for (i in 1:5) {
  for (j in 1:5) {
    set_i <- integration_results[[mechanisms[i]]]
    set_j <- integration_results[[mechanisms[j]]]
    intersection <- sum(set_i & set_j)
    union <- sum(set_i | set_j)
    jaccard_matrix[i, j] <- ifelse(union > 0, intersection / union, 0)
  }
}

pdf(file.path(OUTPUT_DIR, "mechanism_jaccard_heatmap.pdf"), width = 8, height = 7)
pheatmap(
  jaccard_matrix,
  display_numbers = TRUE,
  number_format = "%.3f",
  color = colorRampPalette(c("white", "#2171B5"))(100),
  main = "Mechanism Co-occurrence (Jaccard Index)",
  fontsize = 12,
  fontsize_number = 10
)
dev.off()

cat("\nStep 12: Analyzing genes by mechanism count...\n")

by_n_mech <- integration_results[, .(
  n_genes = .N,
  example_genes = paste(head(gene_symbol[!is.na(gene_symbol)], 5), collapse = ", ")
), by = n_mechanisms][order(n_mechanisms)]

cat("\nGenes by number of mechanisms:\n")
print(by_n_mech)

fwrite(by_n_mech, file.path(OUTPUT_DIR, "genes_by_mechanism_count.csv"))

cat("\nStep 13: Identifying high-priority genes...\n")

if (nrow(multi_mech_genes) > 0) {
  cat("\nTop multi-mechanism genes:\n")
  print(head(multi_mech_genes[, .(gene_symbol, gene_id_clean, n_mechanisms, mechanism_combo)], 20))
}

all_five <- integration_results[n_mechanisms == 5]
if (nrow(all_five) > 0) {
  cat(sprintf("\nGenes with ALL 5 mechanisms: %d\n", nrow(all_five)))
  print(all_five[, .(gene_symbol, gene_id_clean)])
}

four_mech <- integration_results[n_mechanisms == 4]
if (nrow(four_mech) > 0) {
  cat(sprintf("\nGenes with 4 mechanisms: %d\n", nrow(four_mech)))
  print(four_mech[, .(gene_symbol, gene_id_clean, mechanism_combo)])
}

cat("\n============================================\n")
cat("Integration Analysis Complete!\n")
cat("============================================\n\n")

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Output files:\n")
cat("  1. pca_ase_genes_mechanism_integration.csv - Full integration results\n")
cat("  2. mechanism_summary_stats.csv - Summary statistics\n")
cat("  3. mechanism_cooccurrence_matrix.csv - Co-occurrence matrix\n")
cat("  4. mechanism_association_fisher_test.csv - Fisher's exact test results\n")
cat("  5. genes_by_mechanism_count.csv - Genes grouped by mechanism count\n")
cat("  6. genes_*.csv - Genes by individual mechanism\n")
cat("  7. genes_unknown_mechanism.csv - Genes with no identified mechanism\n")
cat("  8. genes_multiple_mechanisms.csv - Genes with >=2 mechanisms\n")
cat("  9. mechanism_barplot.pdf - Bar plot visualization\n")
cat("  10. mechanism_upset.pdf - UpSet plot visualization\n")
cat("  11. primary_mechanism_pie.pdf - Pie chart visualization\n")
cat("  12. mechanism_jaccard_heatmap.pdf - Co-occurrence heatmap\n")

cat("\n============================================\n")
cat("Key Findings:\n")
cat("============================================\n")
cat(sprintf("Total PCa-ASE genes: %d\n", nrow(integration_results)))
cat(sprintf("Genes with identified mechanism: %d (%.1f%%)\n", 
            sum(integration_results$n_mechanisms > 0),
            100 * sum(integration_results$n_mechanisms > 0) / nrow(integration_results)))
cat(sprintf("Genes with multiple mechanisms: %d (%.1f%%)\n",
            sum(integration_results$n_mechanisms >= 2),
            100 * sum(integration_results$n_mechanisms >= 2) / nrow(integration_results)))
cat(sprintf("Genes with unknown mechanism: %d (%.1f%%)\n",
            sum(integration_results$n_mechanisms == 0),
            100 * sum(integration_results$n_mechanisms == 0) / nrow(integration_results)))

cat("\nDone!\n")

cat("\nStep 14: GWAS functional mechanism analysis...\n")

gwas_genes <- integration_results[is_GWAS == TRUE]

if(nrow(gwas_genes) > 0) {
  cat(sprintf("\nGWAS-related genes: %d\n", nrow(gwas_genes)))

  gwas_eqtl <- gwas_genes[is_eQTL == TRUE]
  cat(sprintf("  GWAS + eQTL (expression mechanism): %d (%.1f%%)\n", 
              nrow(gwas_eqtl), 100 * nrow(gwas_eqtl) / nrow(gwas_genes)))

  gwas_sqtl <- gwas_genes[is_sQTL == TRUE]
  cat(sprintf("  GWAS + sQTL (splicing mechanism): %d (%.1f%%)\n",
              nrow(gwas_sqtl), 100 * nrow(gwas_sqtl) / nrow(gwas_genes)))

  gwas_both_qtl <- gwas_genes[is_eQTL == TRUE & is_sQTL == TRUE]
  cat(sprintf("  GWAS + eQTL + sQTL (both mechanisms): %d (%.1f%%)\n",
              nrow(gwas_both_qtl), 100 * nrow(gwas_both_qtl) / nrow(gwas_genes)))

  gwas_only <- gwas_genes[is_eQTL == FALSE & is_sQTL == FALSE]
  cat(sprintf("  GWAS only (unknown mechanism): %d (%.1f%%)\n",
              nrow(gwas_only), 100 * nrow(gwas_only) / nrow(gwas_genes)))

  integration_results[, gwas_mechanism := fifelse(
    !is_GWAS, NA_character_,
    fifelse(is_eQTL & is_sQTL, "GWAS+eQTL+sQTL",
    fifelse(is_eQTL, "GWAS+eQTL",
    fifelse(is_sQTL, "GWAS+sQTL", "GWAS_only")))
  )]

  gwas_mechanism_summary <- integration_results[is_GWAS == TRUE, .N, by = gwas_mechanism]
  fwrite(gwas_mechanism_summary, file.path(OUTPUT_DIR, "gwas_functional_mechanism.csv"))

  gwas_with_qtl <- integration_results[is_GWAS == TRUE & (is_eQTL == TRUE | is_sQTL == TRUE)]
  if(nrow(gwas_with_qtl) > 0) {
    fwrite(gwas_with_qtl, file.path(OUTPUT_DIR, "gwas_genes_with_qtl_support.csv"))
    cat(sprintf("\n  Saved: gwas_genes_with_qtl_support.csv (%d genes)\n", nrow(gwas_with_qtl)))

    cat("\n  High-confidence GWAS genes (with QTL support):\n")
    print(gwas_with_qtl[, .(gene_symbol, gene_id_clean, gwas_mechanism, 
                           has_tumor_eqtl, has_normal_eqtl, 
                           has_tumor_sqtl, has_normal_sqtl)])
  }
}

cat("\nStep 15: SNP-level GWAS-QTL colocalization...\n")

if(file.exists(GWAS_PCA_FILE) && file.exists(QTL_ANNOTATION_FILE)) {

  gwas_snps <- fread(GWAS_PCA_FILE)
  qtl_snps <- fread(QTL_ANNOTATION_FILE)

  if(!"pos_id" %in% names(gwas_snps)) {
    gwas_snps[, pos_id := paste0(gsub("chr", "", seqtnames), ":", start)]
  }

  gwas_sig_snps <- gwas_snps[pvalue < 5e-8, .(pos_id, pvalue)]
  gwas_sig_snps <- unique(gwas_sig_snps, by = "pos_id")

  cat(sprintf("  GWAS significant SNPs (p<5e-8): %d\n", nrow(gwas_sig_snps)))

  qtl_snps_unique <- unique(qtl_snps[, .(
    pos_id, gene_id,
    is_tumor_eqtl, is_normal_eqtl,
    is_tumor_sqtl, is_normal_sqtl
  )])

  gwas_qtl_overlap <- merge(gwas_sig_snps, qtl_snps_unique, by = "pos_id", all.x = TRUE)

  n_gwas_total <- nrow(gwas_sig_snps)
  n_gwas_eqtl <- sum(gwas_qtl_overlap$is_tumor_eqtl | gwas_qtl_overlap$is_normal_eqtl, na.rm = TRUE)
  n_gwas_sqtl <- sum(gwas_qtl_overlap$is_tumor_sqtl | gwas_qtl_overlap$is_normal_sqtl, na.rm = TRUE)

  cat(sprintf("\n  SNP-level colocalization:\n"))
  cat(sprintf("    GWAS SNPs that are also eQTL: %d / %d (%.1f%%)\n",
              n_gwas_eqtl, n_gwas_total, 100 * n_gwas_eqtl / n_gwas_total))
  cat(sprintf("    GWAS SNPs that are also sQTL: %d / %d (%.1f%%)\n",
              n_gwas_sqtl, n_gwas_total, 100 * n_gwas_sqtl / n_gwas_total))

  coloc_snps <- gwas_qtl_overlap[
    (is_tumor_eqtl == TRUE | is_normal_eqtl == TRUE | 
     is_tumor_sqtl == TRUE | is_normal_sqtl == TRUE)
  ]

  if(nrow(coloc_snps) > 0) {
    fwrite(coloc_snps, file.path(OUTPUT_DIR, "gwas_qtl_colocalized_snps.csv"))
    cat(sprintf("    Saved: gwas_qtl_colocalized_snps.csv (%d SNPs)\n", nrow(coloc_snps)))
  }
}
