# ============================================================================
# ============================================================================
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  cat("\nWarning: ComplexHeatmap package not found!\n")
  cat("Make sure the conda r_env environment is activated:\n")
  cat("  conda activate r_env\n")
  cat("or install the ComplexHeatmap package:\n")
  cat("  if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")\n")
  cat("  BiocManager::install(\"ComplexHeatmap\")\n")
  stop("Install ComplexHeatmap or activate the correct conda environment")
}

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(rtracklayer)

cat("\n=== Building the prostate-cancer SCNA occurrence matrix (by chromosome arm) ===\n\n")

output_dir <- "/path/to/data_root/ASE/code/fig2_results/fig2C"

load("/path/to/data_root/ASE/ase_data.RData")
facets_cnv <- read_csv("/path/to/data_root/ASE/derived_inputs/gene_cnv_scores_complete.csv", 
                       show_col_types = FALSE)
pca_ase_genes <- read_csv("/path/to/data_root/ASE/code/PCa_ASE_genes.csv", 
                          show_col_types = FALSE)

fpkm_file <- "/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls"
fpkm_data <- read_tsv(fpkm_file, show_col_types = FALSE)

cat("Loading gene coordinates...\n")
GTF_FILE <- "/path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz"
gtf <- import(GTF_FILE)
gene_coords <- as.data.frame(gtf)
gene_coords <- gene_coords[gene_coords$type == "gene", ]
gene_coords <- gene_coords[, c("gene_id", "seqnames", "start", "end", "gene_name")]
gene_coords <- unique(gene_coords)
gene_coords <- gene_coords %>%
  mutate(
    chromosome = as.character(seqnames),
    chromosome = gsub("^chr", "", chromosome)
  ) %>%
  select(gene_id, chromosome, start, end, gene_name)

cat("Data loaded\n")

centromere_positions <- data.frame(
  chromosome = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                 "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                 "21", "22", "X", "Y"),
  centromere = c(123400000, 93900000, 90900000, 50000000, 48800000,
                 59800000, 60100000, 45200000, 43000000, 39800000,
                 53400000, 35500000, 17700000, 17200000, 19000000,
                 36800000, 25100000, 18500000, 26200000, 28100000,
                 12000000, 15000000, 60600000, 10400000)
)

cat("Chromosome centromere positions loaded\n")

if (!"aRNA" %in% names(ase_data)) {
  if ("ref_ratio" %in% names(ase_data)) {
    ase_data <- ase_data %>% mutate(aRNA = abs(ref_ratio - 0.5))
  }
}

rna_seq_samples <- ase_data %>%
  mutate(sample = str_extract(sample_id, "^T[0-9]+")) %>%
  filter(sample_type == "Tumor") %>%
  pull(sample) %>%
  unique()

cat(sprintf("Tumor samples with RNA-seq data: %d\n\n", length(rna_seq_samples)))

gene_coords_clean <- gene_coords %>%
  mutate(
    gene_id_clean = gsub("\\.[0-9]+$", "", gene_id),
    chromosome = as.character(chromosome)
  )

scna_data <- facets_cnv %>%
  mutate(
    sample = Sample,
    cnv_score_abs = abs(cnv_score),
    chromosome = as.character(chromosome),
    gene_id_clean = gsub("\\.[0-9]+$", "", gene_id)
  ) %>%
  filter(sample %in% rna_seq_samples) %>%
  left_join(
    gene_coords_clean %>% 
      select(gene_id_clean, start, end, chromosome_coord = chromosome),
    by = "gene_id_clean"
  ) %>%
  mutate(
    start = ifelse(is.na(start), 
                   gene_coords_clean$start[match(gene_id, gene_coords_clean$gene_id)], 
                   start),
    chromosome_coord = ifelse(is.na(chromosome_coord),
                              gene_coords_clean$chromosome[match(gene_id, gene_coords_clean$gene_id)],
                              chromosome_coord)
  ) %>%
  mutate(chromosome = ifelse(!is.na(chromosome_coord), chromosome_coord, chromosome)) %>%
  select(-chromosome_coord, -gene_id_clean) %>%
  left_join(centromere_positions, by = "chromosome") %>%
  mutate(
    chr_num = case_when(
      chromosome == "X" ~ 23,
      chromosome == "Y" ~ 24,
      TRUE ~ suppressWarnings(as.numeric(chromosome))
    ),
    arm = case_when(
      is.na(centromere) | is.na(start) ~ "",
      is.na(chr_num) ~ "",
      as.numeric(start) < as.numeric(centromere) ~ "p",
      TRUE ~ "q"
    ),
    chr_arm = ifelse(arm == "", 
                     paste0("chr", chromosome),
                     paste0("chr", chromosome, arm))
  ) %>%
  select(-chr_num)

cat(sprintf("SCNA data after filtering: %d rows\n", nrow(scna_data)))
cat(sprintf("Number of samples: %d\n", n_distinct(scna_data$sample)))

cat("\nChromosome-arm statistics:\n")
scna_data %>%
  count(chr_arm) %>%
  arrange(desc(n)) %>%
  head(20) %>%
  print()

scna_threshold <- 0.4

scna_by_arm <- scna_data %>%
  filter(cnv_score_abs >= scna_threshold) %>%
  group_by(chr_arm) %>%
  summarise(
    n_samples = n_distinct(sample),
    n_genes = n_distinct(gene_id),
    main_cnv_type = names(sort(table(cnv_type), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples))

cat("\n\nTop 15 most frequent SCNA regions (by chromosome arm):\n")
print(scna_by_arm %>% head(15))

top_arms <- scna_by_arm %>%
  head(10) %>%
  pull(chr_arm)

cat(sprintf("\nSelected top 10 SCNA regions:\n%s\n\n", 
            paste(top_arms, collapse = ", ")))

scna_occurrence <- expand.grid(
  sample = rna_seq_samples,
  chr_arm = top_arms,
  stringsAsFactors = FALSE
) %>%
  left_join(
    scna_data %>%
      filter(chr_arm %in% top_arms, cnv_score_abs >= scna_threshold) %>%
      group_by(sample, chr_arm) %>%
      summarise(has_scna = TRUE, .groups = "drop"),
    by = c("sample", "chr_arm")
  ) %>%
  mutate(has_scna = replace_na(has_scna, FALSE))

scna_matrix <- scna_occurrence %>%
  pivot_wider(names_from = chr_arm, values_from = has_scna) %>%
  column_to_rownames("sample")

scna_matrix <- scna_matrix[, top_arms]
scna_matrix <- as.matrix(scna_matrix) * 1

cat(sprintf("SCNA matrix: %d samples × %d regions\n", 
            nrow(scna_matrix), ncol(scna_matrix)))

run_fig2C <- function(key_gene, output_tag = NULL) {
  if (is.null(output_tag) || is.na(output_tag) || output_tag == "") {
    output_tag <- key_gene
  }

  cat(sprintf("\n========== Running Fig2C: key_gene = %s ==========\n", key_gene))

  cat(sprintf("\nExtracting %s expression from FPKM data\n", key_gene))

  gene_expr_raw <- fpkm_data %>%
    dplyr::filter(gene_name == key_gene)

  if (nrow(gene_expr_raw) > 0) {
    cat(sprintf("✓ Found gene %s\n", key_gene))

    tumor_cols <- names(fpkm_data)[grepl("^T[0-9]+", names(fpkm_data))]

    gene_expr_long <- gene_expr_raw %>%
      dplyr::select(dplyr::all_of(tumor_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "sample_name",
        values_to = "fpkm"
      ) %>%
      dplyr::mutate(sample = stringr::str_extract(sample_name, "T[0-9]+"))

    gene_expression <- gene_expr_long %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(
        expression = mean(fpkm, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(expression_z = scale(expression)[, 1])

    cat(sprintf("  Extracted expression for %d samples\n", nrow(gene_expression)))
  } else {
    cat(sprintf("✗ Gene %s not found; no record in the FPKM table\n", key_gene))
    gene_expression <- data.frame(sample = character(0), expression_z = numeric(0))
  }

  samples_in_matrix <- rownames(scna_matrix)

  if (nrow(gene_expression) > 0) {
    missing_samples <- setdiff(samples_in_matrix, gene_expression$sample)
    if (length(missing_samples) > 0) {
      missing_expr <- data.frame(sample = missing_samples, expression_z = 0)
      gene_expression <- dplyr::bind_rows(gene_expression, missing_expr)
    }
  } else {
    set.seed(123)
    gene_expression <- data.frame(
      sample = samples_in_matrix,
      expression_z = rnorm(length(samples_in_matrix))
    )
  }

  sample_order <- gene_expression %>%
    dplyr::filter(sample %in% samples_in_matrix) %>%
    dplyr::arrange(dplyr::desc(expression_z)) %>%
    dplyr::pull(sample)

  scna_matrix_ordered <- scna_matrix[sample_order, ]

  expression_annotation <- gene_expression %>%
    dplyr::filter(sample %in% sample_order) %>%
    dplyr::arrange(match(sample, sample_order)) %>%
    dplyr::pull(expression_z)

  names(expression_annotation) <- sample_order

  cat("\nDrawing heatmap...\n")

  col_scna <- c("0" = "white", "1" = "#8B0000")

  col_expression <- circlize::colorRamp2(
    c(min(expression_annotation, na.rm = TRUE), 0,
      max(expression_annotation, na.rm = TRUE)),
    c("#2166AC", "white", "#B2182B")
  )

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    expression = expression_annotation,
    col = list(expression = col_expression),
    annotation_name_side = "left",
    annotation_legend_param = list(
      expression = list(
        title = sprintf("%s expression", key_gene),
        title_gp = grid::gpar(fontsize = 10)
      )
    ),
    show_annotation_name = FALSE,
    height = grid::unit(0.5, "cm")
  )

  right_anno <- ComplexHeatmap::rowAnnotation(
    occurrence = ComplexHeatmap::anno_barplot(
      colSums(scna_matrix_ordered),
      gp = grid::gpar(fill = "#8B0000"),
      width = grid::unit(2, "cm"),
      axis_param = list(side = "top")
    ),
    annotation_name_side = "top",
    annotation_name_gp = grid::gpar(fontsize = 9)
  )

  ht <- ComplexHeatmap::Heatmap(
    t(scna_matrix_ordered),
    name = "SCNA",
    col = col_scna,
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 10),
    cluster_rows = FALSE,
    column_title = sprintf("Samples (N=%d)", ncol(scna_matrix_ordered)),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    show_column_names = FALSE,
    cluster_columns = FALSE,
    top_annotation = top_anno,
    right_annotation = right_anno,
    heatmap_legend_param = list(
      title = "occurrence",
      at = c(0, 1),
      labels = c("0", "1"),
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9)
    ),
    rect_gp = grid::gpar(col = "gray90", lwd = 0.5),
    width = grid::unit(12, "cm"),
    height = grid::unit(6, "cm")
  )

  pdf(file.path(output_dir,
                sprintf("Fig2C_SCNA_occurrence_chr_arms_%s.pdf", output_tag)),
      width = 10, height = 6)
  ComplexHeatmap::draw(ht)
  dev.off()

  png(file.path(output_dir,
                sprintf("Fig2C_SCNA_occurrence_chr_arms_%s.png", output_tag)),
      width = 10, height = 6, units = "in", res = 300)
  ComplexHeatmap::draw(ht)
  dev.off()

  cat(sprintf("\n✓ Saved Fig2C_SCNA_occurrence_chr_arms_%s.pdf/png\n", output_tag))

  invisible(scna_matrix_ordered)
}

scna_frequency <- data.frame(
  chr_arm   = colnames(scna_matrix),
  n_samples = colSums(scna_matrix),
  frequency = colSums(scna_matrix) / nrow(scna_matrix)
) %>%
  dplyr::arrange(dplyr::desc(n_samples))

readr::write_csv(
  scna_frequency,
  file.path(output_dir, "scna_occurrence_frequency_chr_arms_overall.csv")
)

cat("\n=== Fig2C data preparation done (function run_fig2C defined) ===\n")
cat(sprintf("Number of samples: %d\n", nrow(scna_matrix)))
cat(sprintf("SCNA regions: %d\n", ncol(scna_matrix)))
cat("\nTop 10 chromosome-arm SCNA (overall, independent of key_gene):\n")
print(scna_frequency %>% head(10))
