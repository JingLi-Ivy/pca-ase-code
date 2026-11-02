#!/usr/bin/env Rscript
# ============================================
# Fig2: association of PCa-ASE genes with CNV
# ============================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(stringr)
library(readr)

cat("==========================================\n")
cat("Fig2: PCa-ASE vs CNV association\n")
cat("==========================================\n\n")

PCA_ASE_GENES_FILE <- "/path/to/data_root/ASE/code/PCa_ASE_genes.csv"
FACETS_CNV_FILE <- "/path/to/data_root/variants/batch_43/facets/annotated_cnv_output.tsv"
ASE_DATA_DIR <- "/path/to/data_root/ASE"

OUTPUT_DIR <- "/path/to/data_root/ASE/code/fig2_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Input files:\n")
cat(sprintf("  PCa-ASE genes: %s\n", PCA_ASE_GENES_FILE))
cat(sprintf("  FACETS CNV: %s\n", FACETS_CNV_FILE))
cat(sprintf("  ASE data directory: %s\n", ASE_DATA_DIR))
cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Step 1: Loading PCa-ASE gene list...\n")

pca_ase_genes <- read_csv(PCA_ASE_GENES_FILE, show_col_types = FALSE)
cat(sprintf("  Number of PCa-ASE genes: %d\n", nrow(pca_ase_genes)))

pca_ase_genes <- pca_ase_genes %>%
  mutate(gene_id_clean = gsub("\\..*", "", gene_id))

cat("  Columns: ", paste(colnames(pca_ase_genes), collapse = ", "), "\n")
cat("\n")

cat("Step 2: Loading FACETS CNV data...\n")

facets_cnv <- fread(FACETS_CNV_FILE, sep = "\t")
cat(sprintf("  Total records: %d\n", nrow(facets_cnv)))
cat(sprintf("  Number of samples: %d\n", n_distinct(facets_cnv$sample_id)))
cat(sprintf("  Number of genes: %d\n", n_distinct(facets_cnv$gene_name)))

cat("\nCNV type distribution:\n")
print(table(facets_cnv$cnv_type))
cat("\n")

facets_cnv <- facets_cnv %>%
  mutate(
    tcn_adj = pmax(tcn, 0.1),
    cnv_score = log2(tcn_adj / 2),
    sig_cnv = abs(tcn - 2) >= 1,
    cnv_direction = case_when(
      tcn < 2 ~ "Loss",
      tcn > 2 ~ "Gain",
      TRUE ~ "Neutral"
    )
  )

cat("CNV score distribution:\n")
print(summary(facets_cnv$cnv_score))
cat("\n")

cat("Step 3: Loading ASE data...\n")

ase_files <- list.files(ASE_DATA_DIR, 
                        pattern = "^T[0-9]+_ase_results\\.csv$",
                        full.names = TRUE)

cat(sprintf("  Found %d tumor ASE files\n", length(ase_files)))

ase_data_list <- list()
for (file in ase_files) {
  sample_id <- str_extract(basename(file), "^T[0-9]+")
  df <- fread(file)
  df$Sample <- sample_id
  ase_data_list[[sample_id]] <- df
}

ase_data <- rbindlist(ase_data_list, fill = TRUE)
cat(sprintf("  Total ASE records: %d\n", nrow(ase_data)))
cat(sprintf("  Number of samples: %d\n", n_distinct(ase_data$Sample)))

ase_data <- ase_data %>%
  mutate(gene_id_clean = gsub("\\..*", "", gene_id))

ase_gene_level <- ase_data %>%
  group_by(Sample, gene_id_clean) %>%
  summarise(
    ase_max = max(abs(a), na.rm = TRUE),
    ase_sig = any(fdr < 0.1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.nan(ase_max), !is.infinite(ase_max))

cat(sprintf("  Gene-level ASE records: %d\n\n", nrow(ase_gene_level)))

cat("Step 4: Merging ASE and CNV data...\n")

gene_symbol_to_id <- pca_ase_genes %>%
  dplyr::select(gene_symbol, gene_id_clean) %>%
  distinct()

facets_gene_level <- facets_cnv %>%
  group_by(sample_id, gene_name) %>%
  summarise(
    tcn_avg = mean(tcn, na.rm = TRUE),
    cnv_score = mean(cnv_score, na.rm = TRUE),
    sig_cnv = any(sig_cnv),
    .groups = "drop"
  ) %>%
  rename(Sample = sample_id)

cat(sprintf("  Gene-level CNV records: %d\n", nrow(facets_gene_level)))

ase_cnv_merged <- ase_gene_level %>%
  left_join(gene_symbol_to_id, by = "gene_id_clean") %>%
  inner_join(facets_gene_level, by = c("Sample", "gene_symbol" = "gene_name"))

cat(sprintf("  Records matched by gene_symbol: %d\n\n", nrow(ase_cnv_merged)))

cat("Step 5: Computing ASE-CNV correlations...\n")

pca_ase_gene_ids <- pca_ase_genes$gene_id_clean
pca_ase_gene_symbols <- pca_ase_genes$gene_symbol

ase_cnv_correlations <- ase_cnv_merged %>%
  filter(gene_id_clean %in% pca_ase_gene_ids) %>%
  group_by(gene_id_clean, gene_symbol) %>%
  summarise(
    n_samples = n(),
    rho = if(n() >= 5 && sd(cnv_score, na.rm = TRUE) > 0 && sd(ase_max, na.rm = TRUE) > 0) {
      cor.test(ase_max, cnv_score, method = "spearman", exact = FALSE)$estimate
    } else NA_real_,
    pval = if(n() >= 5 && sd(cnv_score, na.rm = TRUE) > 0 && sd(ase_max, na.rm = TRUE) > 0) {
      cor.test(ase_max, cnv_score, method = "spearman", exact = FALSE)$p.value
    } else NA_real_,
    mean_ase = mean(ase_max, na.rm = TRUE),
    mean_cnv_score = mean(cnv_score, na.rm = TRUE),
    n_scna = sum(sig_cnv, na.rm = TRUE),
    prop_scna = mean(sig_cnv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(rho))

cat(sprintf("  Genes with computed correlations: %d\n", nrow(ase_cnv_correlations)))

ase_cnv_correlations <- ase_cnv_correlations %>%
  mutate(fdr = p.adjust(pval, method = "BH"))

cat("\nClassifying ASE genes...\n")

ase_cnv_correlations <- ase_cnv_correlations %>%
  mutate(
    category = case_when(
      fdr <= 0.1 ~ "ASE, SCNA",
      TRUE ~ "ASE, non-SCNA"
    )
  )

n_ase_scna <- sum(ase_cnv_correlations$category == "ASE, SCNA")
n_ase_non_scna <- sum(ase_cnv_correlations$category == "ASE, non-SCNA")

cat(sprintf("  ASE, SCNA genes: %d\n", n_ase_scna))
cat(sprintf("  ASE, non-SCNA genes: %d\n", n_ase_non_scna))

n_sig_fdr01 <- sum(ase_cnv_correlations$fdr < 0.1, na.rm = TRUE)
n_sig_fdr005 <- sum(ase_cnv_correlations$fdr < 0.05, na.rm = TRUE)
n_sig_p05 <- sum(ase_cnv_correlations$pval < 0.05, na.rm = TRUE)

cat(sprintf("  Significantly correlated genes (FDR < 0.1): %d (%.1f%%)\n", 
            n_sig_fdr01, 100 * n_sig_fdr01 / nrow(ase_cnv_correlations)))
cat(sprintf("  Significantly correlated genes (FDR < 0.05): %d (%.1f%%)\n", 
            n_sig_fdr005, 100 * n_sig_fdr005 / nrow(ase_cnv_correlations)))
cat(sprintf("  Significantly correlated genes (p < 0.05): %d (%.1f%%)\n\n", 
            n_sig_p05, 100 * n_sig_p05 / nrow(ase_cnv_correlations)))

write_csv(ase_cnv_correlations, 
          file.path(OUTPUT_DIR, "pca_ase_cnv_correlations.csv"))

cat("Step 6: Drawing volcano plot...\n")

fdr_threshold_pval <- ase_cnv_correlations %>%
  filter(fdr <= 0.1) %>%
  summarise(max_pval = max(pval, na.rm = TRUE)) %>%
  pull(max_pval)

if (is.na(fdr_threshold_pval) || length(fdr_threshold_pval) == 0) {
  fdr_threshold_pval <- 0.05
}

volcano_data <- ase_cnv_correlations %>%
  mutate(
    neg_log_p = -log10(pval),
    significance = case_when(
      fdr < 0.01 ~ "FDR < 0.01",
      fdr < 0.05 ~ "FDR < 0.05",
      fdr < 0.1 ~ "FDR < 0.1",
      TRUE ~ "n.s."
    ),
    significance = factor(significance, 
                          levels = c("FDR < 0.01", "FDR < 0.05", "FDR < 0.1", "n.s."))
  )

annotation_text <- sprintf("ASE, SCNA: %d\nASE, non-SCNA: %d", 
                           n_ase_scna, n_ase_non_scna)

p_volcano <- ggplot(volcano_data, aes(x = rho, y = neg_log_p, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(fdr_threshold_pval), linetype = "dashed", 
             color = "red", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("FDR < 0.01" = "red", 
               "FDR < 0.05" = "orange", 
               "FDR < 0.1" = "blue", 
               "n.s." = "gray70")
  ) +
  annotate("text", 
           x = min(volcano_data$rho) + 0.05, 
           y = max(volcano_data$neg_log_p) * 0.95,
           label = annotation_text,
           hjust = 0, vjust = 1, size = 4, fontface = "bold") +
  annotate("text", x = max(volcano_data$rho) * 0.8, 
           y = -log10(fdr_threshold_pval) + 0.2,
           label = "FDR = 0.1", color = "red", size = 3, hjust = 0) +
  labs(
    x = "Spearman ρ",
    y = "-log10(P-value)",
    color = "Significance"
  ) +
  theme_classic() +
  theme(
    legend.position = "right"
  )

ggsave(file.path(OUTPUT_DIR, "Fig2A_ASE_CNV_volcano.pdf"), 
       p_volcano, width = 8, height = 6)
ggsave(file.path(OUTPUT_DIR, "Fig2A_ASE_CNV_volcano.png"), 
       p_volcano, width = 8, height = 6, dpi = 300)

cat("  Volcano plot saved\n\n")

cat("Step 6b: Drawing Manhattan plot...\n")

gene_positions <- facets_cnv %>%
  group_by(gene_name) %>%
  summarise(
    chrom = first(chrom),
    start = min(start, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(chrom), chrom %in% c(as.character(1:22), "X", "Y"))

manhattan_data <- ase_cnv_correlations %>%
  left_join(gene_positions, by = c("gene_symbol" = "gene_name")) %>%
  filter(!is.na(chrom)) %>%
  mutate(
    chr_num = case_when(
      chrom == "X" ~ 23,
      chrom == "Y" ~ 24,
      TRUE ~ as.numeric(chrom)
    ),
    neg_log_p = -log10(pval),
    significant = fdr < 0.1
  ) %>%
  arrange(chr_num, start)

cat(sprintf("  Genes with chromosome positions: %d\n", nrow(manhattan_data)))

if (nrow(manhattan_data) > 0) {
  chr_lengths <- manhattan_data %>%
    group_by(chr_num) %>%
    summarise(
      n_genes = n(),
      max_pos = max(start, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(chr_num) %>%
    mutate(
      chr_offset = cumsum(lag(n_genes, default = 0)) + (chr_num - 1) * 5,
      chr_center = chr_offset + n_genes / 2
    )

  manhattan_data <- manhattan_data %>%
    group_by(chr_num) %>%
    mutate(gene_rank = row_number()) %>%
    ungroup() %>%
    left_join(chr_lengths %>% dplyr::select(chr_num, chr_offset), by = "chr_num") %>%
    mutate(x_pos = chr_offset + gene_rank)

  fdr_line <- -log10(fdr_threshold_pval)

  n_scna_manhattan <- sum(manhattan_data$category == "ASE, SCNA")
  n_non_scna_manhattan <- sum(manhattan_data$category == "ASE, non-SCNA")
  annotation_manhattan <- sprintf("ASE, SCNA: %d\nASE, non-SCNA: %d", 
                                  n_scna_manhattan, n_non_scna_manhattan)

  p_manhattan <- ggplot(manhattan_data, aes(x = x_pos, y = neg_log_p)) +
    geom_rect(
      data = chr_lengths %>% filter(chr_num %% 2 == 0),
      aes(xmin = chr_offset, xmax = chr_offset + n_genes, ymin = -Inf, ymax = Inf),
      fill = "gray95", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = fdr_line, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_point(
      data = manhattan_data %>% filter(!significant),
      aes(color = factor(chr_num %% 2)),
      size = 1.5, alpha = 0.6
    ) +
    geom_point(
      data = manhattan_data %>% filter(significant),
      color = "red", size = 2.5, alpha = 0.8
    ) +
    annotate("text", 
             x = max(manhattan_data$x_pos) * 0.02, 
             y = max(manhattan_data$neg_log_p) * 0.95,
             label = annotation_manhattan,
             hjust = 0, vjust = 1, size = 4, fontface = "bold") +
    scale_color_manual(values = c("0" = "gray50", "1" = "gray30"), guide = "none") +
    scale_x_continuous(
      breaks = chr_lengths$chr_center,
      labels = c(1:22, "X", "Y")[1:nrow(chr_lengths)]
    ) +
    labs(
      x = "Chromosome",
      y = "-log10(P-value)"
    ) +
    annotate("text", x = max(manhattan_data$x_pos) * 1.05, y = fdr_line + 0.3,
             label = "FDR = 0.1", color = "red", size = 3, hjust = 1) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 8),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
    )

  ggsave(file.path(OUTPUT_DIR, "Fig2D_Manhattan_plot.pdf"), 
         p_manhattan, width = 14, height = 6)
  ggsave(file.path(OUTPUT_DIR, "Fig2D_Manhattan_plot.png"), 
         p_manhattan, width = 14, height = 6, dpi = 300)

  cat("  Manhattan plot saved\n\n")
} else {
  cat("  Warning: not enough gene-position data for a Manhattan plot\n\n")
}

cat("Step 7: Drawing P-value and rho distributions...\n")

p_pval_hist <- ggplot(ase_cnv_correlations, aes(x = pval)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", 
           x = 0.7, y = Inf, 
           label = sprintf("n = %d genes", nrow(ase_cnv_correlations)),
           hjust = 0, vjust = 1.5, size = 4) +
  labs(
    x = "P-value",
    y = "Count"
  ) +
  theme_classic()

ggsave(file.path(OUTPUT_DIR, "Fig2_pval_distribution.pdf"), 
       p_pval_hist, width = 7, height = 5)

p_rho_hist <- ggplot(ase_cnv_correlations, aes(x = rho)) +
  geom_histogram(bins = 50, fill = "darkgreen", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", 
           x = max(ase_cnv_correlations$rho) * 0.6, y = Inf, 
           label = sprintf("n = %d genes", nrow(ase_cnv_correlations)),
           hjust = 0, vjust = 1.5, size = 4) +
  labs(
    x = "Spearman ρ",
    y = "Count"
  ) +
  theme_classic()

ggsave(file.path(OUTPUT_DIR, "Fig2_rho_distribution.pdf"), 
       p_rho_hist, width = 7, height = 5)

cat("  Distribution plots saved\n\n")

cat("Step 8: Drawing scatter plots of significant genes...\n")

SCNA_SCORE_THRESHOLD <- 0.1

top_genes <- ase_cnv_correlations %>%
  filter(fdr < 0.1) %>%
  arrange(fdr) %>%
  head(6)

if (nrow(top_genes) > 0) {
  cat(sprintf("  Drawing scatter plots for %d significant genes\n", nrow(top_genes)))

  scatter_plots <- list()
  for (i in 1:nrow(top_genes)) {
    gene_id <- top_genes$gene_id_clean[i]
    gene_sym <- top_genes$gene_symbol[i]

    plot_data <- ase_cnv_merged %>%
      filter(gene_id_clean == gene_id) %>%
      mutate(
        point_color = ifelse(ase_sig, "ASE", "no ASE"),
        point_shape = ifelse(abs(cnv_score) >= SCNA_SCORE_THRESHOLD, "SCNA", "CNN")
      )

    if (nrow(plot_data) >= 5) {
      cor_result <- cor.test(plot_data$ase_max, plot_data$cnv_score, 
                             method = "spearman", exact = FALSE)

      gene_fdr <- top_genes$fdr[i]

      p <- ggplot(plot_data, aes(x = ase_max, y = cnv_score, 
                                 color = point_color, shape = point_shape)) +
        geom_point(size = 2.5, alpha = 0.8) +
        geom_hline(yintercept = SCNA_SCORE_THRESHOLD, linetype = "dashed", color = "black") +
        geom_hline(yintercept = -SCNA_SCORE_THRESHOLD, linetype = "dashed", color = "black") +
        scale_color_manual(values = c("ASE" = "red", "no ASE" = "steelblue")) +
        scale_shape_manual(values = c("SCNA" = 17, "CNN" = 16)) +
        annotate(
          "text",
          x = max(plot_data$ase_max, na.rm = TRUE),
          y = max(plot_data$cnv_score, na.rm = TRUE),
          label = sprintf("rho = %.3g\nFDR = %.2g", cor_result$estimate, gene_fdr),
          hjust = 1, vjust = 1, size = 3.5
        ) +
        labs(
          x = expression(ASE(a[RNA])),
          y = "SCNA score",
          color = NULL,
          shape = NULL,
          title = gene_sym
        ) +
        theme_classic(base_size = 10) +
        theme(
          legend.position = "right",
          legend.key.size = unit(0.4, "cm")
        )

      scatter_plots[[i]] <- p

      ggsave(file.path(OUTPUT_DIR, sprintf("scatter_%s.pdf", gene_sym)),
             p, width = 5, height = 4)
    }
  }

  if (length(scatter_plots) >= 2) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      library(patchwork)
      n_plots <- min(6, length(scatter_plots))
      ncol_use <- min(3, n_plots)
      nrow_use <- ceiling(n_plots / ncol_use)
      combined <- wrap_plots(scatter_plots[1:n_plots], 
                             ncol = ncol_use, nrow = nrow_use,
                             guides = "collect")
      ggsave(file.path(OUTPUT_DIR, "Fig2B_top_genes_scatter.pdf"),
             combined, width = 5 * ncol_use, height = 4 * nrow_use)
      ggsave(file.path(OUTPUT_DIR, "Fig2B_top_genes_scatter.png"),
             combined, width = 5 * ncol_use, height = 4 * nrow_use, dpi = 300)
      cat("  Combined scatter plot saved\n")
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      library(gridExtra)
      n_plots <- min(6, length(scatter_plots))
      combined <- grid.arrange(grobs = scatter_plots[1:n_plots],
                               ncol = min(3, n_plots))
      ggsave(file.path(OUTPUT_DIR, "Fig2B_top_genes_scatter.pdf"),
             combined, width = 12, height = 8)
      cat("  Combined scatter plot saved (gridExtra)\n")
    } else {
      cat("  Per-gene scatter plots saved separately\n")
    }
  }
}

cat("\n")

cat("Step 8c: Drawing ATM/INTS6/DSE/MCM9 scatter plots...\n")

target_genes <- c("ATM", "INTS6", "DSE", "MCM9","PVT1","GATA2","NIPAL2","ATP2C1")

for (gene_sym in target_genes) {
  plot_data <- ase_cnv_merged %>%
    filter(gene_symbol == gene_sym) %>%
    mutate(
      ase_abs = ase_max,
      point_color = ifelse(ase_sig, "ASE", "no ASE")
    )

  if (nrow(plot_data) == 0) {
    cat(sprintf(" Warning: no data for %s, skipping\n", gene_sym))
    next
  }

  cor_val <- NA
  cor_p <- NA
  cor_res <- tryCatch(
    cor.test(plot_data$ase_abs, plot_data$cnv_score, method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  if (!is.null(cor_res)) {
    cor_val <- cor_res$estimate
    cor_p <- cor_res$p.value
  }

  gene_fdr <- ase_cnv_correlations %>%
    filter(gene_symbol == gene_sym) %>%
    pull(fdr)
  gene_fdr_txt <- if (length(gene_fdr) > 0) sprintf("FDR = %.2g", gene_fdr[1]) else "FDR = NA"

  p_gene <- ggplot(plot_data, aes(x = ase_abs, y = cnv_score, color = point_color)) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_hline(yintercept = SCNA_SCORE_THRESHOLD, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -SCNA_SCORE_THRESHOLD, linetype = "dashed", color = "black") +
    scale_color_manual(values = c("ASE" = "red", "no ASE" = "blue")) +
    annotate(
      "text",
      x = max(plot_data$ase_abs, na.rm = TRUE),
      y = max(plot_data$cnv_score, na.rm = TRUE),
      label = sprintf("rho = %.3g\n%s", cor_val, gene_fdr_txt),
      hjust = 1, vjust = 1, size = 3.5
    ) +
    labs(
      title = gene_sym,
      x = expression(ASE(a[RNA])),
      y = "SCNA score",
      color = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.4, "cm")
    )

  ggsave(file.path(OUTPUT_DIR, sprintf("scatter_%s_target.pdf", gene_sym)),
         p_gene, width = 5, height = 4)
  ggsave(file.path(OUTPUT_DIR, sprintf("scatter_%s_target.png", gene_sym)),
         p_gene, width = 5, height = 4, dpi = 300)

  cat(sprintf("%s scatter plot saved\n", gene_sym))
}

cat("\n")

cat("Step 8d: Drawing scatter plot of the gene with the largest rho and an SCNA...\n")

high_rho_gene <- ase_cnv_correlations %>%
  filter(prop_scna > 0, rho > 0) %>%
  arrange(desc(rho)) %>%
  slice_head(n = 1)

if (nrow(high_rho_gene) == 0) {
  cat("  Warning: no gene with rho>0 and an SCNA found, skipping\n\n")
} else {
  gene_sym <- high_rho_gene$gene_symbol[1]
  gene_id <- high_rho_gene$gene_id_clean[1]
  gene_rho <- high_rho_gene$rho[1]
  gene_fdr <- high_rho_gene$fdr[1]
  gene_p <- high_rho_gene$pval[1]

  plot_data <- ase_cnv_merged %>%
    filter(gene_id_clean == gene_id) %>%
    mutate(
      point_color = ifelse(ase_sig, "ASE", "no ASE"),
      point_shape = ifelse(abs(cnv_score) >= SCNA_SCORE_THRESHOLD, "SCNA", "CNN")
    )

  if (nrow(plot_data) == 0) {
    cat(sprintf("  Warning: no sample data for %s, skipping\n\n", gene_sym))
  } else {
    p_hr <- ggplot(plot_data, aes(x = ase_max, y = cnv_score,
                                  color = point_color, shape = point_shape)) +
      geom_point(size = 2.5, alpha = 0.85) +
      geom_hline(yintercept = SCNA_SCORE_THRESHOLD, linetype = "dashed", color = "black") +
      geom_hline(yintercept = -SCNA_SCORE_THRESHOLD, linetype = "dashed", color = "black") +
      scale_color_manual(values = c("ASE" = "red", "no ASE" = "steelblue")) +
      scale_shape_manual(values = c("SCNA" = 17, "CNN" = 16)) +
      annotate(
        "text",
        x = max(plot_data$ase_max, na.rm = TRUE),
        y = max(plot_data$cnv_score, na.rm = TRUE),
        label = sprintf("rho = %.3g\np = %.2g\nFDR = %.2g", gene_rho, gene_p, gene_fdr),
        hjust = 1, vjust = 1, size = 3.5
      ) +
      labs(
        title = sprintf("%s (high rho, with SCNA)", gene_sym),
        x = expression(ASE(a[RNA])),
        y = "SCNA score",
        color = NULL,
        shape = NULL
      ) +
      theme_classic(base_size = 10) +
      theme(
        legend.position = "right",
        legend.key.size = unit(0.4, "cm")
      )

    ggsave(file.path(OUTPUT_DIR, sprintf("Fig2B_high_rho_scna_%s.pdf", gene_sym)),
           p_hr, width = 5, height = 4)
    ggsave(file.path(OUTPUT_DIR, sprintf("Fig2B_high_rho_scna_%s.png", gene_sym)),
           p_hr, width = 5, height = 4, dpi = 300)

    cat(sprintf("  Plotted high-rho SCNA gene: %s\n\n", gene_sym))
  }
}

cat("==========================================\n")
cat("Summary\n")
cat("==========================================\n\n")

cat(sprintf("Total PCa-ASE genes: %d\n", nrow(pca_ase_genes)))
cat(sprintf("Genes with ASE-CNV data: %d\n", nrow(ase_cnv_correlations)))
cat(sprintf("  - ASE, SCNA: %d\n", n_ase_scna))
cat(sprintf("  - ASE, non-SCNA: %d\n", n_ase_non_scna))
cat(sprintf("Significant (FDR < 0.1): %d (%.1f%%)\n", 
            n_sig_fdr01, 100 * n_sig_fdr01 / nrow(ase_cnv_correlations)))
cat(sprintf("Positively correlated genes (ρ > 0): %d\n", 
            sum(ase_cnv_correlations$rho > 0, na.rm = TRUE)))
cat(sprintf("Negatively correlated genes (ρ < 0): %d\n", 
            sum(ase_cnv_correlations$rho < 0, na.rm = TRUE)))

cat("\nOutput files:\n")
cat(sprintf("  - %s/pca_ase_cnv_correlations.csv\n", OUTPUT_DIR))
cat(sprintf("  - %s/Fig2A_ASE_CNV_volcano.pdf/png\n", OUTPUT_DIR))
cat(sprintf("  - %s/Fig2B_top_genes_scatter.pdf/png\n", OUTPUT_DIR))
cat(sprintf("  - %s/Fig2D_Manhattan_plot.pdf/png\n", OUTPUT_DIR))
cat(sprintf("  - %s/Fig2_*.pdf (distribution plots)\n", OUTPUT_DIR))

cat("\nAnalysis complete!\n")
