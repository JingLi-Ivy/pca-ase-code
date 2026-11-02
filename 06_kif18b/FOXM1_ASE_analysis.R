#!/usr/bin/env Rscript
# ============================================================================
# FOXM1 expression in samples with KIF18B/NEK2/TOP2A ASE events
# ============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

BASE_DIR <- "/path/to/data_root/ASE/figures"
RNA_FILE <- "/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/Count/DESeq2_normalizedCount_134.xls"
OUTPUT_DIR <- "/path/to/data_root/ASE/figures/FOXM1_analysis"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

target_genes <- c("KIF18B", "NEK2", "TOP2A")

cat(">>> Reading ASE data...\n")

ase_samples <- list()
for (gene in target_genes) {
  snp_file <- file.path(BASE_DIR, paste0(gene, "_analysis"), "All_significant_ASE_SNPs_detail.csv")
  if (file.exists(snp_file)) {
    snp_data <- read_csv(snp_file, show_col_types = FALSE)
    ase_samples[[gene]] <- unique(snp_data$sample_id)
    cat(sprintf("  %s: %d samples with ASE\n", gene, length(ase_samples[[gene]])))
  }
}

samples_with_any_ase <- unique(unlist(ase_samples))
cat(sprintf("\nTotal: %d samples with ASE in any of KIF18B/NEK2/TOP2A\n", length(samples_with_any_ase)))

cat("\n>>> Reading RNA expression data...\n")

rna_raw <- read.delim(RNA_FILE, row.names = 1, check.names = FALSE)
rna_data <- rna_raw %>% tibble::rownames_to_column(var = "gene_name")

foxm1_row <- rna_data %>% filter(gene_name == "FOXM1")

if (nrow(foxm1_row) == 0) {
  stop("FOXM1 gene not found!")
}

foxm1_long <- foxm1_row %>%
  select(-gene_name) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "expression_raw") %>%
  mutate(
    sample_id = str_extract(sample, "^[TN][0-9]+"),
    sample_type = str_extract(sample_id, "^[TN]"),
    expression_raw = as.numeric(expression_raw),
    expression = log2(expression_raw + 1)
  ) %>%
  filter(!is.na(sample_id), sample_type == "T", !is.na(expression)) %>%
  select(sample_id, expression_raw, expression)

cat(sprintf("  FOXM1 expression: %d tumor samples\n", nrow(foxm1_long)))

cat("\n>>> Group comparison (Any gene with ASE vs without ASE)...\n")

foxm1_long <- foxm1_long %>%
  mutate(
    ASE_group = ifelse(sample_id %in% samples_with_any_ase, 
                       "With ASE", "Without ASE"),
    ASE_group = factor(ASE_group, levels = c("Without ASE", "With ASE"))
  )

group_stats <- foxm1_long %>%
  group_by(ASE_group) %>%
  summarise(
    n = n(),
    mean = mean(expression, na.rm = TRUE),
    median = median(expression, na.rm = TRUE),
    sd = sd(expression, na.rm = TRUE),
    .groups = "drop"
  )
print(group_stats)

wilcox_test <- wilcox.test(expression ~ ASE_group, data = foxm1_long)
cat(sprintf("\nWilcoxon test p-value: %.4f\n", wilcox_test$p.value))

p1 <- ggboxplot(foxm1_long, x = "ASE_group", y = "expression",
                fill = "ASE_group", palette = c("#00AFBB", "#FC4E07"),
                add = "jitter", add.params = list(alpha = 0.5, size = 1.5)) +
  stat_compare_means(method = "wilcox.test", label.x = 1.3) +
  labs(
    title = "FOXM1 Expression in Samples with/without\nKIF18B/NEK2/TOP2A ASE Events",
    subtitle = sprintf("With ASE: n=%d; Without ASE: n=%d", 
                       sum(foxm1_long$ASE_group == "With ASE"),
                       sum(foxm1_long$ASE_group == "Without ASE")),
    x = "", y = "FOXM1 Expression (log2)"
  ) +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "FOXM1_vs_any_ASE.pdf"), p1, width = 5, height = 5)

cat("\n>>> Per-gene analysis...\n")

gene_results <- list()

for (gene in target_genes) {
  foxm1_gene <- foxm1_long %>%
    mutate(
      has_ase = ifelse(sample_id %in% ase_samples[[gene]], "With ASE", "Without ASE"),
      has_ase = factor(has_ase, levels = c("Without ASE", "With ASE"))
    )

  stats <- foxm1_gene %>%
    group_by(has_ase) %>%
    summarise(
      n = n(),
      median = median(expression),
      mean = mean(expression),
      .groups = "drop"
    )

  test <- wilcox.test(expression ~ has_ase, data = foxm1_gene)

  gene_results[[gene]] <- list(
    stats = stats,
    p_value = test$p.value
  )

  cat(sprintf("  %s: p=%.4f (With ASE median=%.2f, n=%d; Without ASE median=%.2f, n=%d)\n",
              gene, test$p.value,
              stats$median[stats$has_ase == "With ASE"], stats$n[stats$has_ase == "With ASE"],
              stats$median[stats$has_ase == "Without ASE"], stats$n[stats$has_ase == "Without ASE"]))
}

plot_list <- list()
for (gene in target_genes) {
  foxm1_gene <- foxm1_long %>%
    mutate(
      has_ase = ifelse(sample_id %in% ase_samples[[gene]], "With ASE", "Without ASE"),
      has_ase = factor(has_ase, levels = c("Without ASE", "With ASE"))
    )

  p <- ggboxplot(foxm1_gene, x = "has_ase", y = "expression",
                 fill = "has_ase", palette = c("#00AFBB", "#FC4E07"),
                 add = "jitter", add.params = list(alpha = 0.5)) +
    stat_compare_means(method = "wilcox.test", label.x = 1.2) +
    labs(title = gene, x = "", y = "FOXM1 (log2)") +
    theme(legend.position = "none")

  plot_list[[gene]] <- p
}

p_combined <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 1)
ggsave(file.path(OUTPUT_DIR, "FOXM1_vs_each_gene_ASE.pdf"), p_combined, width = 12, height = 4)

cat("\n>>> FOXM1 high vs low expression groups: KIF18B ASE frequency (CPGEA)...\n")

median_foxm1 <- median(foxm1_long$expression, na.rm = TRUE)

foxm1_groups <- foxm1_long %>%
  mutate(
    FOXM1_group = ifelse(expression >= median_foxm1, "High", "Low"),
    FOXM1_group = factor(FOXM1_group, levels = c("Low", "High")),
    KIF18B_ASE = sample_id %in% ase_samples[["KIF18B"]]
  )

kif18b_freq <- foxm1_groups %>%
  group_by(FOXM1_group) %>%
  summarise(
    n_total = n(),
    n_KIF18B_ASE = sum(KIF18B_ASE),
    freq_KIF18B_ASE = n_KIF18B_ASE / n_total,
    .groups = "drop"
  )
print(kif18b_freq)

write_csv(kif18b_freq,
          file.path(OUTPUT_DIR, "FOXM1_KIF18B_ASE_frequency.csv"))

cont_tab <- table(foxm1_groups$FOXM1_group, foxm1_groups$KIF18B_ASE)
cat("\nContingency table (FOX M1 group vs KIF18B ASE):\n")
print(cont_tab)
fisher_res <- fisher.test(cont_tab)
cat(sprintf("\nFisher's exact test p-value: %.4g\n", fisher_res$p.value))

kif18b_freq_plot <- kif18b_freq %>%
  mutate(
    FOXM1_group = factor(FOXM1_group, levels = c("Low", "High")),
    freq_percent = 100 * freq_KIF18B_ASE
  )

p_kif18b <- ggplot(kif18b_freq_plot,
                   aes(x = FOXM1_group, y = freq_percent, fill = FOXM1_group)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Low" = "#00AFBB", "High" = "#FC4E07"),
                    name = "FOXM1 group") +
  ylim(0, max(kif18b_freq_plot$freq_percent) * 1.25) +
  labs(
    title = "KIF18B ASE frequency by FOXM1 expression",
    subtitle = sprintf("Fisher's exact test p = %.3g", fisher_res$p.value),
    x = "FOXM1 expression group",
    y = "KIF18B ASE frequency (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0, size = 9)
  )

ggsave(
  file.path(OUTPUT_DIR, "KIF18B_ASE_by_FOXM1_group_barplot.pdf"),
  p_kif18b,
  width = 4.5,
  height = 4
)

write_csv(foxm1_long, file.path(OUTPUT_DIR, "FOXM1_expression_with_ASE_status.csv"))

summary_df <- data.frame(
  gene = c("Any", target_genes),
  p_value = c(wilcox_test$p.value, sapply(gene_results, function(x) x$p_value)),
  n_with_ASE = c(sum(foxm1_long$ASE_group == "With ASE"),
                 sapply(target_genes, function(g) length(ase_samples[[g]]))),
  median_with_ASE = c(
    median(foxm1_long$expression[foxm1_long$ASE_group == "With ASE"]),
    sapply(target_genes, function(g) {
      median(foxm1_long$expression[foxm1_long$sample_id %in% ase_samples[[g]]])
    })
  ),
  median_without_ASE = c(
    median(foxm1_long$expression[foxm1_long$ASE_group == "Without ASE"]),
    sapply(target_genes, function(g) {
      median(foxm1_long$expression[!foxm1_long$sample_id %in% ase_samples[[g]]])
    })
  )
)
write_csv(summary_df, file.path(OUTPUT_DIR, "FOXM1_ASE_comparison_summary.csv"))

cat(sprintf("\nResults saved to: %s\n", OUTPUT_DIR))
cat("\nDone!\n")
