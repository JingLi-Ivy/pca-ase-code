#!/usr/bin/env Rscript
# ============================================================
# Unphased methylation extension: methylation-expression-ASE triangulation
#   The base analysis only correlates promoter methylation with ASE (2.8%). This script extends it to a
#   "methylation -> expression -> ASE" triangulation and:
#     (1) includes promoter and gene-body/flanking DMRs (proxy for enhancers)
#     (2) computes methylation-expression and methylation-ASE correlations across tumors
#     (3) flags direction-consistent "methylation-silenced" genes (high methylation -> low expression and -> high ASE)
#     (4) links to the copy_neutral_regulatory events of script 06 to give the
#         "methylation-attributable copy-neutral ASE" subset
#   Note: this is associative (not allele-phased) evidence; label it "methylation-associated"
#         rather than "allele-specific methylation".
#
#   All inputs are read-only; results are written to methylation_triangulation_results/
# ============================================================

suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
})

cat("==========================================\n")
cat("Methylation-expression-ASE triangulation (unphased)\n")
cat("==========================================\n\n")

PCA_ASE_GENES_FILE <- "/path/to/data_root/ASE/code/PCa_ASE_genes.csv"
GTF_FILE   <- "/path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz"
DMR_HYPER  <- "/path/to/data_root/CPGEA/CPGEA_web/WGBS/DMR/table1_hyperDMRs_189pairs_union.txt.gz"
DMR_HYPO   <- "/path/to/data_root/CPGEA/CPGEA_web/WGBS/DMR/table1_hypoDMRs_189pairs_union.txt.gz"
ASE_DIR    <- "/path/to/data_root/ASE"
FPKM_GLOB  <- "/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM*"
EVENT_FILE <- "/path/to/data_root/ASE/code/v8_round2/event_attribution_domseg/event_level_attribution.csv"

OUTPUT_DIR <- "/path/to/data_root/ASE/code/methylation_triangulation_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

PROMOTER_FLANK <- 2000L
MIN_SAMPLES    <- 20L
ASE_FDR        <- 0.10

cat("Step 1: Reading GTF gene coordinates...\n")
gtf <- import(GTF_FILE)
genes <- gtf[gtf$type == "gene"]
gdt <- data.table(
  gene_id = gsub("\\..*", "", genes$gene_id),
  chr = as.character(seqnames(genes)),
  start = start(genes), end = end(genes),
  strand = as.character(strand(genes))
)
gdt <- unique(gdt, by = "gene_id")
gdt[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
gdt[, tss := ifelse(strand == "-", end, start)]
gdt[, prom_start := pmax(1L, tss - PROMOTER_FLANK)]
gdt[, prom_end := tss + PROMOTER_FLANK]
cat(sprintf("  Number of genes: %d\n\n", nrow(gdt)))

pca <- fread(PCA_ASE_GENES_FILE)
pca[, gene_id_clean := gsub("\\..*", "", gene_id)]
gdt <- gdt[gene_id %in% pca$gene_id_clean]
cat(sprintf("  PCa-ASE genes with coordinates: %d\n\n", nrow(gdt)))

cat("Step 2: Loading DMR methylation matrix...\n")
load_dmr <- function(f, type) {
  d <- fread(f)
  setnames(d, 1, "DMR")
  d[, dmr_type := type]
  d
}
dmr <- rbindlist(list(load_dmr(DMR_HYPER, "hyper"),
                      load_dmr(DMR_HYPO,  "hypo")), use.names = TRUE, fill = TRUE)
dmr[, c("dchr","drange") := tstrsplit(DMR, ":", fixed = TRUE)]
dmr[, c("dstart","dend") := tstrsplit(drange, "-", fixed = TRUE)]
dmr[, `:=`(dstart = as.integer(dstart), dend = as.integer(dend))]
dmr[, dchr := ifelse(grepl("^chr", dchr), dchr, paste0("chr", dchr))]
samp_cols <- setdiff(colnames(dmr), c("DMR","dmr_type","dchr","drange","dstart","dend"))
tumor_cols <- grep("^T[0-9]+$", samp_cols, value = TRUE)
cat(sprintf("  DMR regions: %d  tumor sample columns: %d\n\n", nrow(dmr), length(tumor_cols)))

cat("Step 3: Mapping DMRs to gene promoter / gene-body...\n")
dmr[, dmr_idx := .I]
dmr_gr <- GRanges(dmr$dchr, IRanges(dmr$dstart, dmr$dend))
prom_gr <- GRanges(gdt$chr, IRanges(gdt$prom_start, gdt$prom_end), gene_id = gdt$gene_id)
body_gr <- GRanges(gdt$chr, IRanges(gdt$start, gdt$end), gene_id = gdt$gene_id)

map_overlap <- function(query_gr, subj_gr, region_label) {
  hits <- findOverlaps(query_gr, subj_gr)
  if (length(hits) == 0) return(data.table())
  data.table(dmr_idx = queryHits(hits),
             gene_id = subj_gr$gene_id[subjectHits(hits)],
             region = region_label)
}
m_prom <- map_overlap(dmr_gr, prom_gr, "promoter")
m_body <- map_overlap(dmr_gr, body_gr, "genebody")
dmr_gene <- rbindlist(list(m_prom, m_body), fill = TRUE)
dmr_gene[, region := factor(region, levels = c("promoter","genebody"))]
dmr_gene <- dmr_gene[order(dmr_idx, gene_id, region)]
dmr_gene <- unique(dmr_gene, by = c("dmr_idx","gene_id"))
cat(sprintf("  DMR-gene mappings: %d (promoter %d, genebody %d)\n\n",
            nrow(dmr_gene), sum(dmr_gene$region=="promoter"),
            sum(dmr_gene$region=="genebody")))

meth_long <- dmr[dmr_gene, on = "dmr_idx", allow.cartesian = TRUE]
meth_mat <- melt(meth_long, id.vars = c("gene_id","region"),
                 measure.vars = tumor_cols,
                 variable.name = "Sample", value.name = "meth")
meth_gene <- meth_mat[!is.na(meth), .(meth = mean(meth, na.rm = TRUE)),
                      by = .(gene_id, Sample)]
cat(sprintf("  Gene x sample methylation records: %d\n", nrow(meth_gene)))
meth_gene_prom <- meth_mat[region == "promoter" & !is.na(meth),
                           .(meth = mean(meth, na.rm = TRUE)), by = .(gene_id, Sample)]
cat(sprintf("  Of which promoter-only gene x sample: %d\n\n", nrow(meth_gene_prom)))

cat("Step 4: Loading tumor ASE...\n")
ase_files <- list.files(ASE_DIR, pattern = "^T[0-9]+_ase_results\\.csv$", full.names = TRUE)
ase <- rbindlist(lapply(ase_files, function(f) {
  d <- fread(f); d[, Sample := str_extract(basename(f), "^T[0-9]+")]; d
}), fill = TRUE)
ase[, gene_id := gsub("\\..*", "", gene_id)]
ase_gene <- ase[, .(ase_mag = max(abs(a), na.rm = TRUE),
                    ase_sig = any(!is.na(fdr) & fdr < ASE_FDR)),
                by = .(gene_id, Sample)]
ase_gene <- ase_gene[is.finite(ase_mag)]
cat(sprintf("  ASE gene x sample: %d\n\n", nrow(ase_gene)))

cat("Step 5: Loading expression FPKM (expression arm)...\n")
expr_gene <- NULL
fpkm_file <- Sys.glob(FPKM_GLOB)
fpkm_file <- fpkm_file[!grepl("\\.tbi$", fpkm_file)][1]
if (length(fpkm_file) >= 1 && !is.na(fpkm_file) && file.exists(fpkm_file)) {
  tryCatch({
    em <- fread(fpkm_file)
    setnames(em, 1, "gene_id")
    em[, gene_id := gsub("\\..*", "", gene_id)]
    ecols <- grep("^T[0-9]+_WTS$", colnames(em), value = TRUE)
    stopifnot(length(ecols) > 0)
    el <- melt(em, id.vars = "gene_id", measure.vars = ecols,
               variable.name = "Sample", value.name = "fpkm")
    el[, Sample := sub("_WTS$", "", Sample)]
    el[, fpkm := as.numeric(fpkm)]
    expr_gene <- el[is.finite(fpkm)]
    cat(sprintf("  Expression file: %s  tumor sample columns: %d\n", basename(fpkm_file), length(ecols)))
  }, error = function(e) cat("  [warning] failed to load expression:", conditionMessage(e), "\n"))
} else {
  cat("  [warning] FPKM file not found, skipping the expression arm\n")
}
cat("\n")

cat("Step 6: Computing methylation-ASE and methylation-expression correlations per gene...\n")
genes_test <- unique(meth_gene$gene_id)
res <- rbindlist(lapply(genes_test, function(g) {
  mg <- meth_gene[gene_id == g]
  ag <- ase_gene[gene_id == g]
  d <- merge(mg, ag, by = "Sample")
  out <- data.table(gene_id = g, n_meth_ase = nrow(d),
                    rho_meth_ase = NA_real_, p_meth_ase = NA_real_,
                    n_meth_expr = 0L, rho_meth_expr = NA_real_, p_meth_expr = NA_real_,
                    mean_meth = mean(mg$meth, na.rm = TRUE))
  if (nrow(d) >= MIN_SAMPLES && sd(d$meth) > 0 && sd(d$ase_mag) > 0) {
    ct <- suppressWarnings(cor.test(d$meth, d$ase_mag, method = "spearman", exact = FALSE))
    out$rho_meth_ase <- ct$estimate; out$p_meth_ase <- ct$p.value
  }
  dp <- merge(meth_gene_prom[gene_id == g], ag, by = "Sample")
  out$n_prom <- nrow(dp)
  out$rho_meth_ase_prom <- NA_real_; out$p_meth_ase_prom <- NA_real_
  if (nrow(dp) >= MIN_SAMPLES && sd(dp$meth) > 0 && sd(dp$ase_mag) > 0) {
    cp <- suppressWarnings(cor.test(dp$meth, dp$ase_mag, method = "spearman", exact = FALSE))
    out$rho_meth_ase_prom <- cp$estimate; out$p_meth_ase_prom <- cp$p.value
  }
  if (!is.null(expr_gene)) {
    eg <- expr_gene[gene_id == g]
    de <- merge(mg, eg, by = "Sample")
    out$n_meth_expr <- nrow(de)
    if (nrow(de) >= MIN_SAMPLES && sd(de$meth) > 0 && sd(de$fpkm) > 0) {
      ce <- suppressWarnings(cor.test(de$meth, log2(de$fpkm + 1),
                                      method = "spearman", exact = FALSE))
      out$rho_meth_expr <- ce$estimate; out$p_meth_expr <- ce$p.value
    }
  }
  out
}), fill = TRUE)

res[!is.na(p_meth_ase),      fdr_meth_ase      := p.adjust(p_meth_ase, "BH")]
res[!is.na(p_meth_ase_prom), fdr_meth_ase_prom := p.adjust(p_meth_ase_prom, "BH")]
res[!is.na(p_meth_expr),     fdr_meth_expr     := p.adjust(p_meth_expr, "BH")]
res <- merge(res, pca[, .(gene_id = gene_id_clean, gene_symbol)], by = "gene_id", all.x = TRUE)

res[, triangulated := (!is.na(fdr_meth_ase) & fdr_meth_ase < 0.10 & rho_meth_ase > 0 &
                       !is.na(fdr_meth_expr) & fdr_meth_expr < 0.10 & rho_meth_expr < 0)]
fwrite(res[order(-triangulated, fdr_meth_ase)],
       file.path(OUTPUT_DIR, "methylation_triangulation_per_gene.csv"))

n_ma   <- sum(res$fdr_meth_ase < 0.10, na.rm = TRUE)
n_prom <- sum(res$fdr_meth_ase_prom < 0.10, na.rm = TRUE)
n_tri  <- sum(res$triangulated, na.rm = TRUE)
cat(sprintf("\n  Methylation-ASE significant genes (promoter+genebody): %d\n", n_ma))
cat(sprintf("  Of which promoter-only (cf. 66/2.8%% in the base analysis): %d\n", n_prom))
cat(sprintf("  Triangulation-consistent (methylation-silenced) genes: %d\n", n_tri))
if (n_tri > 0) {
  cat("  Triangulation genes (examples):\n")
  print(res[triangulated == TRUE,
            .(gene_symbol, rho_meth_ase, rho_meth_expr, n_meth_expr)][order(rho_meth_expr)][1:min(15,n_tri)])
}

cat("\nStep 7: Linking to copy-neutral events...\n")
attributable_events <- NA_integer_
if (file.exists(EVENT_FILE)) {
  ev <- fread(EVENT_FILE)
  meth_genes <- res[fdr_meth_ase < 0.10 | triangulated == TRUE, gene_id]
  cn <- ev[putative_source == "copy_neutral_regulatory"]
  cn[, meth_attrib := gene_id_clean %in% meth_genes]
  attributable_events <- sum(cn$meth_attrib)
  cat(sprintf("  copy-neutral events: %d, of which methylation-attributable: %d (%.1f%%)\n",
              nrow(cn), attributable_events,
              100*attributable_events/max(1,nrow(cn))))
  fwrite(cn, file.path(OUTPUT_DIR, "copy_neutral_methylation_attributable.csv"))
} else {
  cat("  [note] event-level file not found, skipping (run 06_event_level_attribution.R first)\n")
}

cat("\nStep 8: Plotting...\n")
plot_dt <- res[!is.na(rho_meth_ase) & !is.na(rho_meth_expr)]
if (nrow(plot_dt) > 0) {
  p <- ggplot(plot_dt, aes(x = rho_meth_expr, y = rho_meth_ase, color = triangulated)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_vline(xintercept = 0, color = "grey70") +
    geom_point(size = 1.6, alpha = 0.7) +
    scale_color_manual(values = c(`FALSE`="grey70", `TRUE`="#D6604D"),
                       labels = c("Other","Triangulated"), name = NULL)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = plot_dt[triangulated == TRUE], aes(label = gene_symbol),
      size = 2.8, max.overlaps = 20, show.legend = FALSE)
  }
  p <- p +
    labs(x = "rho(methylation, expression)", y = "rho(methylation, ASE)",
         title = "Methylation-expression-ASE triangulation",
         subtitle = "lower-right = high methylation, low expression, high ASE (methylation-silenced)") +
    theme_bw(base_size = 12)
  ggsave(file.path(OUTPUT_DIR, "Fig_methylation_triangulation.pdf"), p,
         width = 7.5, height = 6)
}

summ <- data.table(
  metric = c("Genes with measurable methylation","Methylation-ASE significant (promoter+genebody, FDR<0.1)",
             "Methylation-ASE significant (promoter-only; base analysis 66/2.8%)","Triangulation-consistent genes",
             "Methylation-attributable copy-neutral events"),
  value = c(length(genes_test), n_ma, n_prom, n_tri, attributable_events)
)
fwrite(summ, file.path(OUTPUT_DIR, "methylation_triangulation_summary.csv"))
cat("\nSummary:\n"); print(summ)

cat("\n==========================================\n")
cat("Done. Output directory:", OUTPUT_DIR, "\n")
cat("  methylation_triangulation_per_gene.csv\n")
cat("  copy_neutral_methylation_attributable.csv\n")
cat("  methylation_triangulation_summary.csv\n")
cat("  Fig_methylation_triangulation.pdf\n")
cat("==========================================\n")
