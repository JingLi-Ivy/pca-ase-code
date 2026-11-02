#!/usr/bin/env Rscript
# ============================================
# PCa-ASE gene mechanism integration, v2
# 
# Framework:
# - ASE mechanisms (4): SCNA, eQTL, sQTL, Methylation
# - GWAS is disease-association evidence, not a mechanism
# - Logic chain: GWAS SNP -> eQTL/sQTL mechanism -> ASE -> disease phenotype
# ============================================

library(data.table)
library(dplyr)
library(ggplot2)
library(UpSetR)

cat("============================================\n")
cat("PCa-ASE Genes Mechanism Integration v2\n")
cat("============================================\n\n")

PCA_ASE_GENES_FILE <- "/path/to/data_root/ASE/code/PCa_ASE_genes.csv"

SCNA_CORR_FILE <- "/path/to/data_root/ASE/code/fig2_results/pca_ase_cnv_correlations.csv"

QTL_ANNOTATION_FILE <- "/path/to/data_root/ASE/code/ASE_SNPs_with_QTL_annotation_exact_match.csv"

GWAS_PCA_FILE <- "/path/to/data_root/ASE/code/PCA_ASE_SNPs_with_GWAS_fixed.csv"

METH_FULL_FILE <- "/path/to/data_root/variants/batch_43/joint_genotyping/pca_ase_methylation_analysis/pca_ase_methylation_analysis_full.txt"

OUTPUT_DIR <- "/path/to/data_root/ASE/code/mechanism_integration_v2"
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

cat("Step 2: Integrating SCNA mechanism...\n")

if (file.exists(SCNA_CORR_FILE)) {
  scna_corr <- fread(SCNA_CORR_FILE)

  if("gene_id_clean" %in% colnames(scna_corr) && "fdr" %in% colnames(scna_corr)) {
    scna_genes <- scna_corr[fdr < 0.1]

    if(nrow(scna_genes) > 0) {
      scna_summary <- scna_genes[, .(
        scna_rho = mean(rho, na.rm = TRUE),
        scna_pval = min(pval, na.rm = TRUE),
        scna_fdr = min(fdr, na.rm = TRUE)
      ), by = gene_id_clean]

      scna_summary[, is_SCNA := TRUE]

      cat(sprintf("  SCNA-driven genes (FDR<0.1): %d\n", nrow(scna_summary)))

      integration_results <- merge(integration_results, scna_summary, 
                                   by = "gene_id_clean", all.x = TRUE)
      integration_results[is.na(is_SCNA), is_SCNA := FALSE]
    } else {
      cat("  No SCNA-driven genes with FDR<0.1\n")
      integration_results[, is_SCNA := FALSE]
    }
  } else {
    cat("  WARNING: Cannot find required columns in SCNA file\n")
    integration_results[, is_SCNA := FALSE]
  }
} else {
  cat("  WARNING: SCNA file not found\n")
  integration_results[, is_SCNA := FALSE]
}

cat("\nStep 3: Integrating eQTL mechanism...\n")

if (file.exists(QTL_ANNOTATION_FILE)) {
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
      eqtl_summary <- eqtl_genes[, .(
        n_eqtl_snps = .N,
        has_tumor_eqtl = if(!is.null(eqtl_tumor_col)) any(get(eqtl_tumor_col) == TRUE, na.rm = TRUE) else FALSE,
        has_normal_eqtl = if(!is.null(eqtl_normal_col)) any(get(eqtl_normal_col) == TRUE, na.rm = TRUE) else FALSE
      ), by = gene_id]

      eqtl_summary[, gene_id_clean := sub("\\..*", "", gene_id)]
      eqtl_summary[, is_eQTL := TRUE]

      cat(sprintf("  eQTL-driven genes: %d\n", nrow(eqtl_summary)))

      integration_results <- merge(
        integration_results, 
        eqtl_summary[, .(gene_id_clean, n_eqtl_snps, has_tumor_eqtl, has_normal_eqtl, is_eQTL)],
        by = "gene_id_clean", all.x = TRUE
      )
      integration_results[is.na(is_eQTL), is_eQTL := FALSE]
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

cat("\nStep 4: Integrating sQTL mechanism...\n")

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
      sqtl_summary <- sqtl_genes[, .(
        n_sqtl_snps = .N,
        has_tumor_sqtl = if(!is.null(sqtl_tumor_col)) any(get(sqtl_tumor_col) == TRUE, na.rm = TRUE) else FALSE,
        has_normal_sqtl = if(!is.null(sqtl_normal_col)) any(get(sqtl_normal_col) == TRUE, na.rm = TRUE) else FALSE
      ), by = gene_id]

      sqtl_summary[, gene_id_clean := sub("\\..*", "", gene_id)]
      sqtl_summary[, is_sQTL := TRUE]

      cat(sprintf("  sQTL-driven genes: %d\n", nrow(sqtl_summary)))

      integration_results <- merge(
        integration_results,
        sqtl_summary[, .(gene_id_clean, n_sqtl_snps, has_tumor_sqtl, has_normal_sqtl, is_sQTL)],
        by = "gene_id_clean", all.x = TRUE
      )
      integration_results[is.na(is_sQTL), is_sQTL := FALSE]
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

cat("\nStep 5: Integrating Methylation mechanism...\n")

if (file.exists(METH_FULL_FILE)) {
  meth_full <- fread(METH_FULL_FILE)

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

      cat(sprintf("  Methylation-driven genes (FDR<0.1): %d\n", nrow(meth_genes)))

      integration_results <- merge(
        integration_results,
        meth_genes,
        by = "gene_id_clean", all.x = TRUE
      )
      integration_results[is.na(is_Methylation), is_Methylation := FALSE]
    } else {
      cat("  No significant methylation-driven genes\n")
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

cat("\nStep 6: Adding GWAS disease association evidence...\n")

integration_results[, `:=`(
  has_GWAS_evidence = FALSE,
  n_gwas_snps = NA_integer_,
  min_gwas_pval = NA_real_
)]

if (file.exists(GWAS_PCA_FILE) && file.exists(QTL_ANNOTATION_FILE)) {

  gwas_data <- fread(GWAS_PCA_FILE)
  qtl_data <- fread(QTL_ANNOTATION_FILE)

  pval_col <- if("pvalue" %in% colnames(gwas_data)) "pvalue" else if("pval" %in% colnames(gwas_data)) "pval" else if("p.value" %in% colnames(gwas_data)) "p.value" else NULL

  if(!is.null(pval_col)) {
    gwas_sig <- gwas_data[get(pval_col) < 5e-8]

    if(nrow(gwas_sig) > 0) {
      cat(sprintf("  GWAS significant SNPs (p<5e-8): %d\n", nrow(gwas_sig)))

      if(!"pos_id" %in% names(gwas_sig)) {
        gwas_sig[, pos_id := paste0(gsub("chr", "", seqnames), ":", start)]
      }

      if("pos_id" %in% names(qtl_data) && "gene_id" %in% names(qtl_data)) {
        qtl_snp_map <- unique(qtl_data[!is.na(gene_id), .(pos_id, gene_id)])
        qtl_snp_map[, gene_id_clean := sub("\\..*", "", gene_id)]

        gwas_genes <- merge(
          gwas_sig[, .(pos_id, gwas_pval = get(pval_col))],
          qtl_snp_map[, .(pos_id, gene_id_clean)],
          by = "pos_id"
        )

        if(nrow(gwas_genes) > 0) {
          gwas_summary <- gwas_genes[, .(
            n_gwas_snps = .N,
            min_gwas_pval = min(gwas_pval, na.rm = TRUE)
          ), by = gene_id_clean]

          gwas_summary[, has_GWAS_evidence := TRUE]

          cat(sprintf("  Genes with GWAS evidence: %d\n", nrow(gwas_summary)))

          integration_results <- merge(
            integration_results[, !c("has_GWAS_evidence", "n_gwas_snps", "min_gwas_pval"), with = FALSE],
            gwas_summary,
            by = "gene_id_clean", all.x = TRUE
          )
          integration_results[is.na(has_GWAS_evidence), has_GWAS_evidence := FALSE]
        }
      }
    } else {
      cat("  No GWAS significant SNPs found\n")
    }
  }
} else {
  cat("  WARNING: GWAS or QTL file not found\n")
}

cat("\nStep 7: Summarizing mechanisms...\n")

integration_results[, n_mechanisms := 
                      as.integer(is_SCNA) + 
                      as.integer(is_eQTL) + 
                      as.integer(is_sQTL) + 
                      as.integer(is_Methylation)
]

integration_results[, mechanism_combo := paste0(
  ifelse(is_SCNA, "SCNA,", ""),
  ifelse(is_eQTL, "eQTL,", ""),
  ifelse(is_sQTL, "sQTL,", ""),
  ifelse(is_Methylation, "Meth,", "")
)]
integration_results[, mechanism_combo := sub(",$", "", mechanism_combo)]
integration_results[mechanism_combo == "", mechanism_combo := "Unknown"]

integration_results[, primary_mechanism := fifelse(
  is_SCNA, "SCNA",
  fifelse(is_eQTL, "eQTL",
          fifelse(is_sQTL, "sQTL",
                  fifelse(is_Methylation, "Methylation", "Unknown"))))
]

integration_results[, gwas_functional_class := fifelse(
  !has_GWAS_evidence, "No_GWAS",
  fifelse(is_eQTL & is_sQTL, "GWAS_eQTL_sQTL",
          fifelse(is_eQTL, "GWAS_eQTL",
                  fifelse(is_sQTL, "GWAS_sQTL", "GWAS_mechanism_unknown")))
)]

cat("\n============================================\n")
cat("Summary Statistics\n")
cat("============================================\n\n")

cat("Total PCa-ASE genes:", nrow(integration_results), "\n\n")

cat("[ASE Mechanisms (4 types)]\n")
cat(sprintf("  SCNA-driven:        %d (%.1f%%)\n", 
            sum(integration_results$is_SCNA), 
            100 * sum(integration_results$is_SCNA) / nrow(integration_results)))
cat(sprintf("  eQTL-driven:        %d (%.1f%%)\n",
            sum(integration_results$is_eQTL),
            100 * sum(integration_results$is_eQTL) / nrow(integration_results)))
cat(sprintf("  sQTL-driven:        %d (%.1f%%)\n",
            sum(integration_results$is_sQTL),
            100 * sum(integration_results$is_sQTL) / nrow(integration_results)))
cat(sprintf("  Methylation-driven: %d (%.1f%%)\n",
            sum(integration_results$is_Methylation),
            100 * sum(integration_results$is_Methylation) / nrow(integration_results)))
cat(sprintf("  Unknown mechanism:  %d (%.1f%%)\n",
            sum(integration_results$n_mechanisms == 0),
            100 * sum(integration_results$n_mechanisms == 0) / nrow(integration_results)))

cat("\n[GWAS Disease Association Evidence]\n")
cat(sprintf("  Genes with GWAS evidence: %d (%.1f%%)\n",
            sum(integration_results$has_GWAS_evidence),
            100 * sum(integration_results$has_GWAS_evidence) / nrow(integration_results)))

cat("\n[High-confidence causal genes (GWAS + QTL mechanism)]\n")
gwas_eqtl <- sum(integration_results$has_GWAS_evidence & integration_results$is_eQTL)
gwas_sqtl <- sum(integration_results$has_GWAS_evidence & integration_results$is_sQTL)
gwas_both <- sum(integration_results$has_GWAS_evidence & integration_results$is_eQTL & integration_results$is_sQTL)
gwas_unknown <- sum(integration_results$has_GWAS_evidence & !integration_results$is_eQTL & !integration_results$is_sQTL)

cat(sprintf("  GWAS + eQTL (expression mechanism): %d\n", gwas_eqtl))
cat(sprintf("  GWAS + sQTL (splicing mechanism):   %d\n", gwas_sqtl))
cat(sprintf("  GWAS + eQTL + sQTL (both):          %d\n", gwas_both))
cat(sprintf("  GWAS only (mechanism unclear):      %d\n", gwas_unknown))

cat("\nGenes by number of mechanisms:\n")
print(table(integration_results$n_mechanisms))

cat("\nStep 8: Generating visualizations...\n")

mechanism_counts <- data.frame(
  Mechanism = c("SCNA", "eQTL", "sQTL", "Methylation", "Unknown"),
  Count = c(
    sum(integration_results$is_SCNA),
    sum(integration_results$is_eQTL),
    sum(integration_results$is_sQTL),
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
    x = "ASE Mechanism",
    y = "Number of Genes"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  ylim(0, max(mechanism_counts$Count) * 1.15)

ggsave(file.path(OUTPUT_DIR, "mechanism_barplot.pdf"), p1, width = 10, height = 6)

upset_data <- integration_results[, .(
  SCNA = as.integer(is_SCNA),
  eQTL = as.integer(is_eQTL),
  sQTL = as.integer(is_sQTL),
  Methylation = as.integer(is_Methylation)
)]

pdf(file.path(OUTPUT_DIR, "mechanism_upset.pdf"), width = 10, height = 7)
upset(
  as.data.frame(upset_data),
  sets = c("SCNA", "eQTL", "sQTL", "Methylation"),
  order.by = "freq",
  decreasing = TRUE,
  mb.ratio = c(0.6, 0.4),
  text.scale = 1.5,
  main.bar.color = "#4292C6",
  sets.bar.color = "#2171B5",
  matrix.color = "#2171B5"
)
dev.off()

gwas_func_counts <- integration_results[has_GWAS_evidence == TRUE, .N, by = gwas_functional_class]
gwas_func_counts[, pct := 100 * N / sum(N)]

if(nrow(gwas_func_counts) > 0) {
  p3 <- ggplot(gwas_func_counts, aes(x = reorder(gwas_functional_class, -N), y = N, fill = gwas_functional_class)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%d (%.1f%%)", N, pct)), vjust = -0.3, size = 3.5) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "GWAS Genes by Functional Mechanism",
      subtitle = "Causal mechanism for GWAS-associated ASE genes",
      x = "Functional Class",
      y = "Number of Genes"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    ylim(0, max(gwas_func_counts$N) * 1.15)

  ggsave(file.path(OUTPUT_DIR, "gwas_functional_mechanism.pdf"), p3, width = 8, height = 6)
}

cat("  Visualizations saved\n")

cat("\nStep 9: Saving mechanism-specific tables...\n")

base_cols <- c("gene_id", "gene_id_clean", "gene_symbol", 
               "n_mechanisms", "mechanism_combo", "primary_mechanism",
               "has_GWAS_evidence", "n_gwas_snps", "min_gwas_pval", "gwas_functional_class")

scna_genes <- integration_results[is_SCNA == TRUE]
if(nrow(scna_genes) > 0) {
  scna_cols <- c(base_cols, "scna_rho", "scna_pval", "scna_fdr")
  scna_cols <- scna_cols[scna_cols %in% names(scna_genes)]
  scna_output <- scna_genes[, ..scna_cols]
  setorder(scna_output, -has_GWAS_evidence, scna_fdr)
  fwrite(scna_output, file.path(OUTPUT_DIR, "sheet1_SCNA_driven_genes.csv"))
  cat(sprintf("  Saved: sheet1_SCNA_driven_genes.csv (%d genes)\n", nrow(scna_output)))
}

eqtl_genes <- integration_results[is_eQTL == TRUE]
if(nrow(eqtl_genes) > 0) {
  eqtl_cols <- c(base_cols, "n_eqtl_snps", "has_tumor_eqtl", "has_normal_eqtl")
  eqtl_cols <- eqtl_cols[eqtl_cols %in% names(eqtl_genes)]
  eqtl_output <- eqtl_genes[, ..eqtl_cols]
  setorder(eqtl_output, -has_GWAS_evidence, -n_eqtl_snps)
  fwrite(eqtl_output, file.path(OUTPUT_DIR, "sheet2_eQTL_driven_genes.csv"))
  cat(sprintf("  Saved: sheet2_eQTL_driven_genes.csv (%d genes)\n", nrow(eqtl_output)))
}

sqtl_genes <- integration_results[is_sQTL == TRUE]
if(nrow(sqtl_genes) > 0) {
  sqtl_cols <- c(base_cols, "n_sqtl_snps", "has_tumor_sqtl", "has_normal_sqtl")
  sqtl_cols <- sqtl_cols[sqtl_cols %in% names(sqtl_genes)]
  sqtl_output <- sqtl_genes[, ..sqtl_cols]
  setorder(sqtl_output, -has_GWAS_evidence, -n_sqtl_snps)
  fwrite(sqtl_output, file.path(OUTPUT_DIR, "sheet3_sQTL_driven_genes.csv"))
  cat(sprintf("  Saved: sheet3_sQTL_driven_genes.csv (%d genes)\n", nrow(sqtl_output)))
}

meth_genes <- integration_results[is_Methylation == TRUE]
if(nrow(meth_genes) > 0) {
  meth_cols <- c(base_cols, "meth_rho", "meth_fdr", "n_dmr", "dmr_types")
  meth_cols <- meth_cols[meth_cols %in% names(meth_genes)]
  meth_output <- meth_genes[, ..meth_cols]
  setorder(meth_output, -has_GWAS_evidence, meth_fdr)
  fwrite(meth_output, file.path(OUTPUT_DIR, "sheet4_Methylation_driven_genes.csv"))
  cat(sprintf("  Saved: sheet4_Methylation_driven_genes.csv (%d genes)\n", nrow(meth_output)))
}

unknown_genes <- integration_results[n_mechanisms == 0]
if(nrow(unknown_genes) > 0) {
  unknown_cols <- base_cols[base_cols %in% names(unknown_genes)]
  unknown_output <- unknown_genes[, ..unknown_cols]
  setorder(unknown_output, -has_GWAS_evidence, gene_symbol)
  fwrite(unknown_output, file.path(OUTPUT_DIR, "sheet5_Unknown_mechanism_genes.csv"))
  cat(sprintf("  Saved: sheet5_Unknown_mechanism_genes.csv (%d genes)\n", nrow(unknown_output)))
}

high_conf_genes <- integration_results[has_GWAS_evidence == TRUE & (is_eQTL == TRUE | is_sQTL == TRUE)]
if(nrow(high_conf_genes) > 0) {
  all_cols <- names(high_conf_genes)
  high_conf_output <- high_conf_genes[, ..all_cols]
  setorder(high_conf_output, min_gwas_pval)
  fwrite(high_conf_output, file.path(OUTPUT_DIR, "sheet6_High_confidence_causal_genes.csv"))
  cat(sprintf("  Saved: sheet6_High_confidence_causal_genes.csv (%d genes)\n", nrow(high_conf_output)))
}

setorder(integration_results, -n_mechanisms, primary_mechanism, gene_symbol)
fwrite(integration_results, file.path(OUTPUT_DIR, "full_integration_results.csv"))
cat(sprintf("  Saved: full_integration_results.csv (%d genes)\n", nrow(integration_results)))

cat("\nStep 10: Saving summary statistics...\n")

summary_stats <- data.table(
  Category = c(
    "Total PCa-ASE genes",
    "--- ASE Mechanisms ---",
    "SCNA-driven",
    "eQTL-driven",
    "sQTL-driven",
    "Methylation-driven",
    "Unknown mechanism",
    "--- Mechanism counts ---",
    "Single mechanism",
    "Multiple mechanisms (>=2)",
    "--- GWAS Evidence ---",
    "Has GWAS association",
    "GWAS + eQTL",
    "GWAS + sQTL",
    "GWAS + eQTL + sQTL",
    "GWAS only (mechanism unclear)"
  ),
  Count = c(
    nrow(integration_results),
    NA,
    sum(integration_results$is_SCNA),
    sum(integration_results$is_eQTL),
    sum(integration_results$is_sQTL),
    sum(integration_results$is_Methylation),
    sum(integration_results$n_mechanisms == 0),
    NA,
    sum(integration_results$n_mechanisms == 1),
    sum(integration_results$n_mechanisms >= 2),
    NA,
    sum(integration_results$has_GWAS_evidence),
    gwas_eqtl,
    gwas_sqtl,
    gwas_both,
    gwas_unknown
  )
)

summary_stats[!is.na(Count), Percentage := sprintf("%.1f%%", 100 * Count / nrow(integration_results))]

fwrite(summary_stats, file.path(OUTPUT_DIR, "summary_statistics.csv"))
cat("  Saved: summary_statistics.csv\n")

cat("\nStep 11: Mechanism co-occurrence analysis...\n")

mechanisms <- c("is_SCNA", "is_eQTL", "is_sQTL", "is_Methylation")
mech_names <- c("SCNA", "eQTL", "sQTL", "Methylation")

cross_tab <- matrix(0, nrow = 4, ncol = 4)
rownames(cross_tab) <- mech_names
colnames(cross_tab) <- mech_names

for (i in 1:4) {
  for (j in 1:4) {
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

if(length(fisher_results) > 0) {
  fisher_dt <- rbindlist(fisher_results)
  fisher_dt[, FDR := p.adjust(pvalue, method = "BH")]
  fisher_dt <- fisher_dt[order(pvalue)]

  cat("\nSignificant associations (FDR < 0.05):\n")
  print(fisher_dt[FDR < 0.05])

  fwrite(fisher_dt, file.path(OUTPUT_DIR, "mechanism_association_fisher_test.csv"))
}

library(pheatmap)

jaccard_matrix <- matrix(0, nrow = 4, ncol = 4)
rownames(jaccard_matrix) <- mech_names
colnames(jaccard_matrix) <- mech_names

for (i in 1:4) {
  for (j in 1:4) {
    set_i <- integration_results[[mechanisms[i]]]
    set_j <- integration_results[[mechanisms[j]]]
    intersection <- sum(set_i & set_j)
    union <- sum(set_i | set_j)
    jaccard_matrix[i, j] <- ifelse(union > 0, intersection / union, 0)
  }
}

pdf(file.path(OUTPUT_DIR, "mechanism_jaccard_heatmap.pdf"), width = 7, height = 6)
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

cat("\n============================================\n")
cat("Integration Analysis Complete!\n")
cat("============================================\n\n")

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Output files for Excel sheets:\n")
cat("  sheet1_SCNA_driven_genes.csv\n")
cat("  sheet2_eQTL_driven_genes.csv\n")
cat("  sheet3_sQTL_driven_genes.csv\n")
cat("  sheet4_Methylation_driven_genes.csv\n")
cat("  sheet5_Unknown_mechanism_genes.csv\n")
cat("  sheet6_High_confidence_causal_genes.csv\n")
cat("\nSupplementary files:\n")
cat("  full_integration_results.csv\n")
cat("  summary_statistics.csv\n")
cat("  mechanism_cooccurrence_matrix.csv\n")
cat("  mechanism_association_fisher_test.csv\n")
cat("\nVisualizations:\n")
cat("  mechanism_barplot.pdf\n")
cat("  mechanism_upset.pdf\n")
cat("  gwas_functional_mechanism.pdf\n")
cat("  mechanism_jaccard_heatmap.pdf\n")

cat("\n============================================\n")
cat("Key Findings:\n")
cat("============================================\n")
cat(sprintf("Total PCa-ASE genes: %d\n", nrow(integration_results)))
cat(sprintf("Genes with identified mechanism: %d (%.1f%%)\n", 
            sum(integration_results$n_mechanisms > 0),
            100 * sum(integration_results$n_mechanisms > 0) / nrow(integration_results)))
cat(sprintf("High-confidence causal genes (GWAS+QTL): %d\n", nrow(high_conf_genes)))

cat("\nDone!\n")
