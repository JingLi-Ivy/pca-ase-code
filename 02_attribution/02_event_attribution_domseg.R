#!/usr/bin/env Rscript
# ============================================================
# Event-level (sample-gene) source attribution + confidence tiers
#   Upgrades non-exclusive gene-level overlap to a per tumor-gene ASE event
#   source attribution (deletion / amplification / copy-neutral LOH /
#   balanced copy-neutral / unresolved) with a confidence tier.
#
#   Inputs (read-only):
#     - PCa-ASE gene list (gene_id <-> gene_symbol map)
#     - per-tumor gene-level ASE results (T*_ase_results.csv)
#     - FACETS-annotated CNV (annotated_cnv_output.tsv)
#     - [optional] non-SCNA source annotations (eQTL/sQTL/methylation/missense)
#   Outputs (written to a new directory; raw data are not modified):
#     - event_level_attribution.csv   one row per ASE event
#     - event_level_source_proportions.csv  event-level source proportions
#     - gene_vs_event_proportions.csv gene-level vs event-level comparison
#     - Fig_event_level_attribution.pdf
# ============================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(readr)
  library(ggplot2)
})

cat("==========================================\n")
cat("Event-level source attribution (sample-gene level)\n")
cat("==========================================\n\n")

PCA_ASE_GENES_FILE <- "/path/to/data_root/ASE/code/PCa_ASE_genes.csv"
FACETS_CNV_FILE    <- "/path/to/data_root/variants/batch_43/facets/annotated_cnv_output.tsv"
ASE_DATA_DIR       <- "/path/to/data_root/ASE"
SOURCE_ANNOT_FILE  <- "/path/to/data_root/ASE/code/Gene_level_ASE_QTL_annotation.csv"

OUTPUT_DIR <- "/path/to/data_root/ASE/code/v8_round2/event_attribution_domseg"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

ASE_FDR    <- 0.10
DEL_TCN    <- 1.5
AMP_TCN    <- 2.5

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Step 1: Loading PCa-ASE gene list...\n")
pca <- read_csv(PCA_ASE_GENES_FILE, show_col_types = FALSE) %>%
  mutate(gene_id_clean = gsub("\\..*", "", gene_id))
gene_map <- pca %>% dplyr::select(gene_id_clean, gene_symbol) %>% distinct()
cat(sprintf("  Number of PCa-ASE genes: %d\n\n", nrow(pca)))

cat("Step 2: Loading tumor ASE results and extracting significant events...\n")
ase_files <- list.files(ASE_DATA_DIR,
                        pattern = "^T[0-9]+_ase_results\\.csv$",
                        full.names = TRUE)
cat(sprintf("  Tumor ASE files: %d\n", length(ase_files)))
stopifnot(length(ase_files) > 0)

ase_list <- lapply(ase_files, function(f) {
  s  <- str_extract(basename(f), "^T[0-9]+")
  d  <- fread(f)
  d[, Sample := s]
  d
})
ase <- rbindlist(ase_list, fill = TRUE)
ase[, gene_id_clean := gsub("\\..*", "", gene_id)]

ase <- ase[gene_id_clean %in% pca$gene_id_clean]
ase[, ase_sig := (!is.na(fdr) & fdr < ASE_FDR)]
ase_evt <- ase[, .(ase_a = max(abs(a), na.rm = TRUE),
                   ase_fdr = min(fdr, na.rm = TRUE),
                   ase_sig = any(ase_sig)),
               by = .(Sample, gene_id_clean)]
ase_evt <- ase_evt[is.finite(ase_a)]
sig_evt <- ase_evt[ase_sig == TRUE]
cat(sprintf("  Significant ASE events (tumor x gene, fdr<%.2f): %d\n", ASE_FDR, nrow(sig_evt)))
cat(sprintf("  Genes involved: %d  samples involved: %d\n\n",
            uniqueN(sig_evt$gene_id_clean), uniqueN(sig_evt$Sample)))

cat("Step 3: Loading FACETS and aggregating tcn/lcn by event key...\n")
fac <- fread(FACETS_CNV_FILE, sep = "\t")

evt_base <- as.data.table(
  sig_evt %>% left_join(gene_map, by = "gene_id_clean")
)
evt_keys <- unique(evt_base[!is.na(gene_symbol), .(sample_id = Sample, gene_name = gene_symbol)])

fac_sub <- fac[evt_keys, on = .(sample_id, gene_name), nomatch = 0]
fac_sub[, tcn := suppressWarnings(as.numeric(tcn))][, lcn := suppressWarnings(as.numeric(lcn))]
setorder(fac_sub, sample_id, gene_name, -num_mark)
fac_gene <- fac_sub[, .(tcn = tcn[1L], lcn = lcn[1L], cnv_type = as.character(cnv_type[1L])),
                    by = .(sample_id, gene_name)]
setnames(fac_gene, c("sample_id", "gene_name"), c("Sample", "gene_symbol"))
cat(sprintf("  FACETS (sample,gene) records for the events: %d\n\n", nrow(fac_gene)))

cat("Step 4: Attributing each ASE event...\n")
evt <- evt_base %>%
  as_tibble() %>%
  left_join(fac_gene, by = c("Sample", "gene_symbol"))

evt <- evt %>%
  mutate(event_class = case_when(
    is.na(tcn)                 ~ "UNRESOLVED",
    tcn < DEL_TCN              ~ "DELETION",
    tcn > AMP_TCN & !is.na(lcn) & abs(2 * lcn - tcn) < 1e-9 ~ "BALANCED_GAIN",
    tcn > AMP_TCN              ~ "AMPLIFICATION",
    !is.na(lcn) & lcn == 0     ~ "COPY_NEUTRAL_LOH",
    !is.na(lcn) & lcn >= 1     ~ "BALANCED",
    TRUE                       ~ "UNRESOLVED"
  ))

evt <- evt %>%
  mutate(putative_source = case_when(
    event_class == "DELETION"          ~ "SCNA_deletion",
    event_class == "AMPLIFICATION"     ~ "SCNA_amplification",
    event_class == "COPY_NEUTRAL_LOH"  ~ "CN_LOH",
    event_class == "BALANCED"          ~ "copy_neutral_regulatory",
    event_class == "BALANCED_GAIN"     ~ "balanced_gain",
    TRUE                               ~ "unresolved"
  ))

has_annot <- FALSE
if (file.exists(SOURCE_ANNOT_FILE)) {
  annot <- tryCatch(fread(SOURCE_ANNOT_FILE), error = function(e) NULL)
  if (!is.null(annot)) {
    nm <- tolower(colnames(annot))
    id_col <- colnames(annot)[which(nm %in% c("gene_id","gene_id_clean","ensembl","gene"))][1]
    if (!is.na(id_col)) {
      annot[, gene_id_clean := gsub("\\..*", "", get(id_col))]
      flag_cols <- colnames(annot)[grepl("eqtl|sqtl|methyl|missense|mutation", nm)]
      if (length(flag_cols) > 0) {
        annot[, nonscna_evidence := rowSums(sapply(.SD, function(x) {
          as.integer(!is.na(x) & x != 0 & x != "" & tolower(as.character(x)) != "false")
        })) > 0, .SDcols = flag_cols]
        evt <- evt %>%
          left_join(annot[, .(gene_id_clean, nonscna_evidence)] %>% distinct(),
                    by = "gene_id_clean")
        has_annot <- TRUE
        cat(sprintf("  Non-SCNA annotation added; source columns: %s\n",
                    paste(flag_cols, collapse = ", ")))
      }
    }
  }
}
if (!has_annot) {
  cat("  [note] no usable non-SCNA source annotation columns found; copy-neutral events are treated as unannotated\n")
  evt$nonscna_evidence <- NA
}

evt <- evt %>%
  mutate(source_confidence = case_when(
    putative_source %in% c("SCNA_deletion","SCNA_amplification","CN_LOH") ~ "HIGH",
    putative_source %in% c("copy_neutral_regulatory", "balanced_gain") &
      !is.na(nonscna_evidence) & nonscna_evidence                         ~ "MODERATE",
    TRUE                                                                  ~ "LOW"
  ))

evt_out <- evt %>%
  dplyr::select(Sample, gene_id_clean, gene_symbol, ase_a, ase_fdr,
                tcn, lcn, cnv_type, event_class, putative_source, source_confidence)
write_csv(evt_out, file.path(OUTPUT_DIR, "event_level_attribution.csv"))
cat(sprintf("  Wrote event-level attribution table: %d rows\n\n", nrow(evt_out)))

cat("Step 5: Summarizing proportions...\n")
evt_prop <- evt_out %>%
  count(putative_source, name = "n_events") %>%
  mutate(pct_events = round(100 * n_events / sum(n_events), 1)) %>%
  arrange(desc(n_events))
write_csv(evt_prop, file.path(OUTPUT_DIR, "event_level_source_proportions.csv"))
cat("Event-level source proportions:\n"); print(evt_prop)

conf_prop <- evt_out %>%
  count(source_confidence, name = "n_events") %>%
  mutate(pct = round(100 * n_events / sum(n_events), 1))
cat("\nConfidence tiers:\n"); print(conf_prop)

gene_level <- evt_out %>%
  group_by(gene_id_clean) %>%
  summarise(
    any_scna = any(putative_source %in% c("SCNA_deletion","SCNA_amplification","CN_LOH")),
    any_cn   = any(putative_source == "copy_neutral_regulatory"),
    .groups = "drop"
  )
gene_vs_event <- bind_rows(
  evt_prop %>% transmute(level = "event", source = putative_source, pct = pct_events),
  gene_level %>%
    summarise(SCNA_any = 100*mean(any_scna),
              copy_neutral_any = 100*mean(any_cn)) %>%
    pivot_longer(everything(), names_to = "source", values_to = "pct") %>%
    transmute(level = "gene", source, pct = round(pct,1))
)
write_csv(gene_vs_event, file.path(OUTPUT_DIR, "gene_vs_event_proportions.csv"))
cat("\nGene-level vs event-level comparison written\n\n")

cat("Step 6: Plotting...\n")
src_levels <- c("SCNA_deletion","SCNA_amplification","CN_LOH",
                "copy_neutral_regulatory","unresolved")
src_cols <- c(SCNA_deletion="#2166AC", SCNA_amplification="#B2182B",
              CN_LOH="#4393C3", copy_neutral_regulatory="#F4A582",
              unresolved="grey70")

p1 <- evt_prop %>%
  mutate(putative_source = factor(putative_source, levels = src_levels)) %>%
  ggplot(aes(x = reorder(putative_source, -n_events), y = n_events, fill = putative_source)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n_events, pct_events)),
            vjust = -0.2, size = 3) +
  scale_fill_manual(values = src_cols, guide = "none") +
  labs(x = NULL, y = "No. of ASE events",
       title = "Event-level source attribution") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p2 <- conf_prop %>%
  mutate(source_confidence = factor(source_confidence,
                                    levels = c("HIGH","MODERATE","LOW"))) %>%
  ggplot(aes(x = source_confidence, y = n_events, fill = source_confidence)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.2, size = 3.2) +
  scale_fill_manual(values = c(HIGH="#1A9850", MODERATE="#FEE08B", LOW="grey75"),
                    guide = "none") +
  labs(x = "Source confidence", y = "No. of ASE events",
       title = "Event-level confidence tier") +
  theme_bw(base_size = 12)

pdf(file.path(OUTPUT_DIR, "Fig_event_level_attribution.pdf"), width = 11, height = 4.5)
if (requireNamespace("gridExtra", quietly = TRUE)) {
  gridExtra::grid.arrange(p1, p2, ncol = 2)
} else {
  print(p1); print(p2)
}
dev.off()

cat("\n==========================================\n")
cat("Done. Output files:\n")
cat(sprintf("  %s/event_level_attribution.csv\n", OUTPUT_DIR))
cat(sprintf("  %s/event_level_source_proportions.csv\n", OUTPUT_DIR))
cat(sprintf("  %s/gene_vs_event_proportions.csv\n", OUTPUT_DIR))
cat(sprintf("  %s/Fig_event_level_attribution.pdf\n", OUTPUT_DIR))
cat("==========================================\n")
