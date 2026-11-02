#!/usr/bin/env Rscript
# Targeted cis-trans interaction screen for recurrent PCa-ASE cis variants.
#
# Candidate cis variants:
#   1. Recurrent significant ASE SNPs in PCa-ASE genes (n_sig_samples >= 5)
#   2. Existing strong LD coding proxies (r2 >= 0.8) for PCa-ASE SNPs
#   3. Forced KIF18B proxy rs9914422, to anchor the validated locus
#
# Model per tissue, candidate variant-gene pair, and regulator:
#   log2(target_gene_FPKM + 0.01) ~ dosage + regulator_z + dosage:regulator_z

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(ggplot2)
})

PROJECT_DIR <- "/path/to/data_root/ASE"
OUTDIR <- Sys.getenv("OUTDIR")
TRANSFORM <- Sys.getenv("TRANSFORM", "INT")
KEEP_FILE <- Sys.getenv("KEEP_FILE", "")
FPKM_FILE <- "/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls"
PLINK_PREFIX <- "/path/to/data_root/CPGEA/CPGEA_genotype/combined.134.qc"
PCA_ASE_GENES_FILE <- file.path(PROJECT_DIR, "code/PCa_ASE_genes.csv")
ASE_FILE <- file.path(PROJECT_DIR, "all_T_samples_SNPs_with_ASE_annotation.csv")
LD_PROXY_FILE <- file.path(PROJECT_DIR, "code/fig4_results/all_2378_pca_ase_snps_vep_with_ldproxy.csv")

MIN_SIG_SAMPLES <- 5L
MIN_N <- 80L
MIN_MAC <- 10
MIN_TARGET_MEDIAN <- 0.1
MIN_REGULATOR_MEDIAN <- 0.1
CHUNK_SIZE <- 80L
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- 1L

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

tf_panel <- unique(c(
  "AR", "FOXA1", "FOXA2", "FOXA3", "HOXB13", "NKX3-1", "GATA2", "GATA3",
  "ERG", "ETV1", "ETV4", "ETV5", "FLI1", "SPI1", "SPDEF", "EHF",
  "ELF1", "ELF3", "ELF4", "ETS1", "ETS2",
  "FOXM1", "FOXK1", "FOXK2", "FOXO1", "FOXO3", "FOXO4", "FOXP1", "FOXP2",
  "FOXP4", "FOXJ1", "FOXC1", "FOXC2", "FOXD1",
  "MYC", "MYCN", "MAX", "MXD1", "MNT", "E2F1", "E2F2", "E2F3", "E2F4",
  "E2F5", "E2F6", "E2F7", "E2F8", "TFDP1", "TFDP2", "MYBL1", "MYBL2",
  "MYB", "TP53", "TP63", "TP73",
  "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "ATF1", "ATF2",
  "ATF3", "ATF4", "ATF5", "ATF6", "ATF7", "CREB1", "CREB3", "CREM",
  "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", "CEBPZ",
  "RELA", "RELB", "NFKB1", "NFKB2", "STAT1", "STAT2", "STAT3", "STAT4",
  "STAT5A", "STAT5B", "STAT6", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5",
  "IRF7", "IRF8", "IRF9",
  "ESR1", "ESR2", "ESRRA", "ESRRB", "NR3C1", "NR3C2", "NR2F1", "NR2F2",
  "PPARA", "PPARG", "RXRA", "RXRB", "RARA", "RARB", "RARG", "VDR",
  "HNF1A", "HNF1B", "HNF4A", "AHR",
  "SOX2", "SOX4", "SOX9", "SOX10", "POU2F1", "POU2F2", "POU2F3",
  "POU5F1", "KLF2", "KLF4", "KLF5", "KLF6", "KLF9", "KLF10", "KLF13",
  "SP1", "SP2", "SP3", "SP4", "GRHL1", "GRHL2", "GRHL3", "RUNX1",
  "RUNX2", "RUNX3", "ZBTB16", "BCL6", "PRDM1", "REST",
  "SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2", "ID1", "ID2",
  "ID3", "ID4", "TCF3", "TCF4", "TCF7", "TCF7L1", "TCF7L2", "TCF12",
  "LEF1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "SMAD2", "SMAD3", "SMAD4",
  "HIF1A", "EPAS1", "HES1", "HES2", "HES5", "HEY1", "HEY2",
  "MEF2A", "MEF2B", "MEF2C", "MEF2D", "SRF", "YAP1", "TAZ",
  "HOXA9", "HOXA10", "HOXA13", "HOXB9", "HOXB13", "HOXC6", "HOXC8",
  "HOXC10", "HOXD10", "HOXD13", "NKX2-1", "NKX2-5",
  "CTCF", "BRD4", "EZH2", "SUZ12", "EED", "CHD1", "CHD4", "SMARCA4",
  "SMARCA2", "ARID1A", "ARID1B", "EP300", "CREBBP"
))

decode_plink_snp <- function(raw_bytes, n_samples) {
  byte_val <- as.integer(raw_bytes)
  bits0 <- byte_val %% 4L
  bits1 <- (byte_val %/% 4L) %% 4L
  bits2 <- (byte_val %/% 16L) %% 4L
  bits3 <- (byte_val %/% 64L) %% 4L
  bits <- as.integer(rbind(bits0, bits1, bits2, bits3))
  bits <- bits[seq_len(n_samples)]
  dosage <- rep(NA_real_, n_samples)
  dosage[bits == 0L] <- 2
  dosage[bits == 2L] <- 1
  dosage[bits == 3L] <- 0
  dosage
}

read_plink_dosage <- function(plink_prefix, variant_indices) {
  fam <- fread(paste0(plink_prefix, ".fam"), header = FALSE)
  setnames(fam, c("fid", "iid", "father", "mother", "sex", "pheno"))
  n_samples <- nrow(fam)
  bytes_per_snp <- ceiling(n_samples / 4)
  variant_indices <- sort(unique(as.integer(variant_indices)))

  con <- file(paste0(plink_prefix, ".bed"), "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, "raw", n = 3)
  if (!identical(magic, as.raw(c(0x6c, 0x1b, 0x01)))) {
    stop("PLINK BED file is not in SNP-major mode or has an unexpected header.")
  }

  mat <- matrix(NA_real_, nrow = n_samples, ncol = length(variant_indices))
  colnames(mat) <- as.character(variant_indices)
  for (j in seq_along(variant_indices)) {
    idx <- variant_indices[j]
    seek(con, where = 3 + (idx - 1L) * bytes_per_snp, origin = "start")
    raw_bytes <- readBin(con, "raw", n = bytes_per_snp)
    mat[, j] <- decode_plink_snp(raw_bytes, n_samples)
  }

  sample_id <- sub("^N", "", fam$fid)
  list(sample_id = sample_id, dosage = mat)
}

collapse_expression_rows <- function(rows, sample_cols) {
  if (nrow(rows) == 0) return(rows)
  rows[, row_median := apply(as.matrix(.SD), 1, median, na.rm = TRUE),
       .SDcols = sample_cols]
  setorder(rows, expr_key, -row_median)
  rows[, .SD[1], by = expr_key]
}

make_expression_wide <- function(fpkm, target_gene_ids, regulator_symbols, tissue_prefix) {
  sample_cols <- grep(paste0("^", tissue_prefix, "\\d+_WTS$"), names(fpkm), value = TRUE)
  target_rows <- fpkm[gene_id_clean %in% target_gene_ids]
  target_rows[, expr_key := gene_id_clean]
  regulator_rows <- fpkm[gene_name %in% regulator_symbols]
  regulator_rows[, expr_key := gene_name]
  rows <- rbindlist(list(target_rows, regulator_rows), use.names = TRUE, fill = TRUE)
  rows <- collapse_expression_rows(rows, sample_cols)

  expr <- melt(
    rows[, c("expr_key", sample_cols), with = FALSE],
    id.vars = "expr_key",
    variable.name = "sample_col",
    value.name = "FPKM"
  )
  expr[, sample_id := sub(paste0("^", tissue_prefix), "", sample_col)]
  expr[, sample_id := sub("_WTS$", "", sample_id)]
  dcast(expr, sample_id ~ expr_key, value.var = "FPKM")
}

int_transform <- function(v) {
  qnorm((rank(v, ties.method = "average") - 0.5) / length(v))
}

fit_interaction_fast <- function(y, x, z) {
  ok <- is.finite(y) & is.finite(x) & is.finite(z)
  n <- sum(ok)
  if (n < MIN_N) return(NULL)
  y0 <- y[ok]
  x0 <- x[ok]
  z0 <- z[ok]
  if (length(unique(x0)) < 2L) return(NULL)
  mac <- min(sum(x0, na.rm = TRUE), 2 * n - sum(x0, na.rm = TRUE))
  if (mac < MIN_MAC) return(NULL)
  if (length(unique(z0)) < 3L) return(NULL)
  if (TRANSFORM == "INT") {
    y0 <- int_transform(y0)
    z0 <- int_transform(z0)
  } else {
    z0 <- as.numeric(scale(z0))
  }
  sst <- sum((y0 - mean(y0))^2)
  if (!is.finite(sst) || sst <= 0) return(NULL)

  xz <- x0 * z0
  x_add <- cbind(1, x0, z0)
  x_int <- cbind(1, x0, z0, xz)
  fit_add <- lm.fit(x_add, y0)
  fit_int <- lm.fit(x_int, y0)
  if (fit_int$rank < 4L || fit_add$rank < 3L) return(NULL)

  rss_add <- sum(fit_add$residuals^2)
  rss_int <- sum(fit_int$residuals^2)
  df_num <- fit_int$rank - fit_add$rank
  df_den <- n - fit_int$rank
  if (df_num < 1L || df_den < 1L || rss_int <= 0) return(NULL)
  f_stat <- max(0, (rss_add - rss_int) / df_num) / (rss_int / df_den)
  p_int <- pf(f_stat, df_num, df_den, lower.tail = FALSE)

  list(
    n = n,
    mac = mac,
    beta_interact = fit_int$coefficients[4],
    p_interact = p_int,
    R2_additive = 1 - rss_add / sst,
    R2_interaction = 1 - rss_int / sst,
    delta_R2 = (rss_add - rss_int) / sst
  )
}

build_candidates <- function() {
  cat("Loading PCa-ASE genes and significant ASE SNPs...\n")
  pca <- fread(PCA_ASE_GENES_FILE)
  pca[, gene_id_clean := sub("\\..*$", "", gene_id)]

  ase <- fread(
    ASE_FILE,
    select = c("seqtnames", "start", "end", "ref", "alt", "snp_id",
               "gene_id", "sample_id", "fdr", "ase_status")
  )
  ase[, gene_id_clean := sub("\\..*$", "", gene_id)]
  ase_sig <- ase[
    gene_id_clean %in% pca$gene_id_clean &
      (ase_status == "Significant_ASE" | (!is.na(fdr) & fdr < 0.1))
  ]
  direct <- ase_sig[, .(
    n_sig_samples = uniqueN(sample_id),
    n_sig_records = .N,
    min_fdr = min(fdr, na.rm = TRUE)
  ), by = .(gene_id_clean, seqtnames, start, end, ref, alt, snp_id)]
  direct <- direct[n_sig_samples >= MIN_SIG_SAMPLES]
  direct <- merge(
    direct,
    unique(pca[, .(gene_id_clean, gene_symbol)]),
    by = "gene_id_clean",
    all.x = TRUE
  )
  direct[, `:=`(
    candidate_source = "direct_recurrent_ase_snp",
    variant_rsid = NA_character_,
    proxy_r2 = NA_real_,
    proxy_func = NA_character_,
    source_snp_id = snp_id
  )]

  cat("Loading existing LD proxy summary...\n")
  ld <- fread(
    LD_PROXY_FILE,
    select = c("gene_symbol", "chr", "pos", "query_rsid", "best_proxy_rsid",
               "best_proxy_r2", "best_proxy_func")
  )
  ld_proxy <- ld[
    !is.na(best_proxy_rsid) & best_proxy_rsid != "" &
      !is.na(best_proxy_r2) & best_proxy_r2 >= 0.8
  ]
  setorder(ld_proxy, gene_symbol, -best_proxy_r2)
  ld_proxy <- ld_proxy[, .SD[1], by = .(gene_symbol, best_proxy_rsid)]
  ld_proxy <- merge(
    ld_proxy,
    unique(pca[, .(gene_id_clean, gene_symbol)]),
    by = "gene_symbol",
    all.x = FALSE
  )
  ld_proxy[, `:=`(
    candidate_source = "ld_coding_proxy_r2_ge_0.8",
    seqtnames = NA_character_,
    start = NA_integer_,
    end = NA_integer_,
    ref = NA_character_,
    alt = NA_character_,
    snp_id = NA_character_,
    variant_rsid = best_proxy_rsid,
    proxy_r2 = best_proxy_r2,
    proxy_func = best_proxy_func,
    source_snp_id = query_rsid,
    n_sig_samples = NA_integer_,
    n_sig_records = NA_integer_,
    min_fdr = NA_real_
  )]

  forced <- data.table(
    gene_id_clean = "ENSG00000186185",
    gene_symbol = "KIF18B",
    seqtnames = NA_character_,
    start = NA_integer_,
    end = NA_integer_,
    ref = NA_character_,
    alt = NA_character_,
    snp_id = NA_character_,
    n_sig_samples = NA_integer_,
    n_sig_records = NA_integer_,
    min_fdr = NA_real_,
    candidate_source = "forced_kif18b_validated_proxy",
    variant_rsid = "rs9914422",
    proxy_r2 = NA_real_,
    proxy_func = "validated_KIF18B_cis_proxy",
    source_snp_id = "rs9914422"
  )

  candidates <- rbindlist(list(
    direct[, names(forced), with = FALSE],
    ld_proxy[, names(forced), with = FALSE],
    forced
  ), use.names = TRUE, fill = TRUE)
  unique(candidates)
}

map_candidates_to_bim <- function(candidates, bim) {
  bim[, bed_index := .I]
  bim[, chr := as.character(chr)]
  bim[, key_a1a2 := paste(chr, pos, a1, a2, sep = "|")]
  bim[, key_a2a1 := paste(chr, pos, a2, a1, sep = "|")]

  direct <- candidates[candidate_source == "direct_recurrent_ase_snp"]
  direct[, chr := as.character(seqtnames)]
  direct[, pos := as.integer(start)]
  direct[, coord_key := paste(chr, pos, ref, alt, sep = "|")]

  m1 <- merge(
    direct,
    bim[, .(key_a1a2, bed_index, bed_snp = snp, bed_chr = chr, bed_pos = pos,
            bed_a1 = a1, bed_a2 = a2)],
    by.x = "coord_key",
    by.y = "key_a1a2",
    all = FALSE
  )
  if (nrow(m1)) {
    m1[, `:=`(
      invert_dosage = TRUE,
      dosage_allele = alt,
      mapping_method = "coord_ref_a1_alt_a2"
    )]
  }

  m2 <- merge(
    direct,
    bim[, .(key_a2a1, bed_index, bed_snp = snp, bed_chr = chr, bed_pos = pos,
            bed_a1 = a1, bed_a2 = a2)],
    by.x = "coord_key",
    by.y = "key_a2a1",
    all = FALSE
  )
  if (nrow(m2)) {
    m2[, `:=`(
      invert_dosage = FALSE,
      dosage_allele = alt,
      mapping_method = "coord_ref_a2_alt_a1"
    )]
  }

  rsid <- candidates[!is.na(variant_rsid) & variant_rsid != ""]
  mr <- merge(
    rsid,
    bim[, .(bed_snp = snp, bed_index, bed_chr = chr, bed_pos = pos,
            bed_a1 = a1, bed_a2 = a2)],
    by.x = "variant_rsid",
    by.y = "bed_snp",
    all = FALSE
  )
  if (nrow(mr)) {
    mr[, bed_snp := variant_rsid]
    mr[, `:=`(
      invert_dosage = FALSE,
      dosage_allele = bed_a1,
      mapping_method = "rsid_a1"
    )]
    mr[variant_rsid == "rs9914422" & bed_a1 == "C" & bed_a2 == "G",
       `:=`(
         invert_dosage = TRUE,
         dosage_allele = "G",
         mapping_method = "rsid_forced_G_dosage"
       )]
  }

  mapped <- rbindlist(list(m1, m2, mr), use.names = TRUE, fill = TRUE)
  mapped[, source_priority := fifelse(
    candidate_source == "forced_kif18b_validated_proxy", 0L,
    fifelse(candidate_source == "direct_recurrent_ase_snp", 1L, 2L)
  )]
  mapped[, recurrence_for_sort := fifelse(is.na(n_sig_samples), -1L, n_sig_samples)]
  setorder(mapped, gene_id_clean, bed_index, source_priority, -recurrence_for_sort)
  mapped <- mapped[, .SD[1], by = .(gene_id_clean, bed_index)]
  mapped[, model_variant_id := fifelse(!is.na(variant_rsid) & variant_rsid != "",
                                       variant_rsid, bed_snp)]
  mapped[, locus_id := paste(gene_id_clean, model_variant_id, sep = "|")]
  mapped
}

prepare_tissue_data <- function(expr_wide, geno, mapped, tissue_label) {
  expr <- expr_wide[match(geno$sample_id, expr_wide$sample_id)]
  sample_ok <- !is.na(expr$sample_id)
  if (nzchar(KEEP_FILE)) sample_ok <- sample_ok & geno$sample_id %in% readLines(KEEP_FILE)
  expr <- expr[sample_ok]
  dosage <- geno$dosage[sample_ok, , drop = FALSE]
  list(expr = expr, dosage = dosage, tissue = tissue_label)
}

run_tissue_screen <- function(tissue_data, mapped, tf_symbols) {
  expr <- tissue_data$expr
  dosage_mat <- tissue_data$dosage
  tissue_label <- tissue_data$tissue

  available_tf <- intersect(tf_symbols, names(expr))
  tf_raw <- as.matrix(expr[, available_tf, with = FALSE])
  tf_median <- apply(tf_raw, 2, median, na.rm = TRUE)
  tf_sd <- apply(tf_raw, 2, sd, na.rm = TRUE)
  keep_tf <- is.finite(tf_median) & tf_median >= MIN_REGULATOR_MEDIAN &
    is.finite(tf_sd) & tf_sd > 0
  available_tf <- available_tf[keep_tf]
  tf_raw <- tf_raw[, available_tf, drop = FALSE]
  tf_median <- tf_median[available_tf]

  target_gene_ids <- intersect(unique(mapped$gene_id_clean), names(expr))
  target_median <- sapply(target_gene_ids, function(g) median(expr[[g]], na.rm = TRUE))
  target_gene_ids <- target_gene_ids[is.finite(target_median) & target_median >= MIN_TARGET_MEDIAN]
  mapped_tissue <- mapped[gene_id_clean %in% target_gene_ids]

  cat(sprintf("[%s] regulators tested: %d; candidate loci with expressed targets: %d\n",
              tissue_label, length(available_tf), nrow(mapped_tissue)))
  if (nrow(mapped_tissue) == 0 || length(available_tf) == 0) return(data.table())

  chunks <- split(seq_len(nrow(mapped_tissue)),
                  ceiling(seq_len(nrow(mapped_tissue)) / CHUNK_SIZE))
  worker <- function(idx) {
    max_rows <- length(idx) * length(available_tf)
    out_i <- 0L
    tissue_v <- character(max_rows)
    gene_id_v <- character(max_rows)
    gene_symbol_v <- character(max_rows)
    model_variant_id_v <- character(max_rows)
    bed_snp_v <- character(max_rows)
    bed_chr_v <- character(max_rows)
    bed_pos_v <- integer(max_rows)
    dosage_allele_v <- character(max_rows)
    candidate_source_v <- character(max_rows)
    source_snp_id_v <- character(max_rows)
    n_sig_samples_v <- integer(max_rows)
    proxy_r2_v <- numeric(max_rows)
    tf_v <- character(max_rows)
    n_v <- integer(max_rows)
    mac_v <- numeric(max_rows)
    n_dosage0_v <- integer(max_rows)
    n_dosage1_v <- integer(max_rows)
    n_dosage2_v <- integer(max_rows)
    median_tf_v <- numeric(max_rows)
    beta_interact_v <- numeric(max_rows)
    p_interact_v <- numeric(max_rows)
    r2_additive_v <- numeric(max_rows)
    r2_interaction_v <- numeric(max_rows)
    delta_r2_v <- numeric(max_rows)
    locus_id_v <- character(max_rows)

    for (ii in idx) {
      cand <- mapped_tissue[ii]
      y <- log2(expr[[cand$gene_id_clean]] + 0.01)
      x <- dosage_mat[, as.character(cand$bed_index)]
      if (isTRUE(cand$invert_dosage)) x <- 2 - x
      g_counts <- table(factor(x, levels = c(0, 1, 2)), useNA = "no")

      for (tf in available_tf) {
        fit <- fit_interaction_fast(y, x, tf_raw[, tf])
        if (is.null(fit)) next
        out_i <- out_i + 1L
        tissue_v[out_i] <- tissue_label
        gene_id_v[out_i] <- cand$gene_id_clean
        gene_symbol_v[out_i] <- cand$gene_symbol
        model_variant_id_v[out_i] <- cand$model_variant_id
        bed_snp_v[out_i] <- cand$bed_snp
        bed_chr_v[out_i] <- cand$bed_chr
        bed_pos_v[out_i] <- cand$bed_pos
        dosage_allele_v[out_i] <- cand$dosage_allele
        candidate_source_v[out_i] <- cand$candidate_source
        source_snp_id_v[out_i] <- cand$source_snp_id
        n_sig_samples_v[out_i] <- ifelse(is.na(cand$n_sig_samples), NA_integer_, cand$n_sig_samples)
        proxy_r2_v[out_i] <- ifelse(is.na(cand$proxy_r2), NA_real_, cand$proxy_r2)
        tf_v[out_i] <- tf
        n_v[out_i] <- fit$n
        mac_v[out_i] <- fit$mac
        n_dosage0_v[out_i] <- as.integer(g_counts[1])
        n_dosage1_v[out_i] <- as.integer(g_counts[2])
        n_dosage2_v[out_i] <- as.integer(g_counts[3])
        median_tf_v[out_i] <- tf_median[tf]
        beta_interact_v[out_i] <- fit$beta_interact
        p_interact_v[out_i] <- fit$p_interact
        r2_additive_v[out_i] <- fit$R2_additive
        r2_interaction_v[out_i] <- fit$R2_interaction
        delta_r2_v[out_i] <- fit$delta_R2
        locus_id_v[out_i] <- cand$locus_id
      }
    }
    if (out_i == 0L) return(data.table())
    keep <- seq_len(out_i)
    data.table(
      tissue = tissue_v[keep],
      gene_id = gene_id_v[keep],
      gene_symbol = gene_symbol_v[keep],
      model_variant_id = model_variant_id_v[keep],
      bed_snp = bed_snp_v[keep],
      bed_chr = bed_chr_v[keep],
      bed_pos = bed_pos_v[keep],
      dosage_allele = dosage_allele_v[keep],
      candidate_source = candidate_source_v[keep],
      source_snp_id = source_snp_id_v[keep],
      n_sig_samples = n_sig_samples_v[keep],
      proxy_r2 = proxy_r2_v[keep],
      TF = tf_v[keep],
      n = n_v[keep],
      mac = mac_v[keep],
      n_dosage0 = n_dosage0_v[keep],
      n_dosage1 = n_dosage1_v[keep],
      n_dosage2 = n_dosage2_v[keep],
      median_TF = median_tf_v[keep],
      beta_interact = beta_interact_v[keep],
      p_interact = p_interact_v[keep],
      R2_additive = r2_additive_v[keep],
      R2_interaction = r2_interaction_v[keep],
      delta_R2 = delta_r2_v[keep],
      locus_id = locus_id_v[keep]
    )
  }

  n_cores <- min(N_CORES, length(chunks))
  if (n_cores <= 1L) {
    results <- rbindlist(lapply(chunks, worker), use.names = TRUE, fill = TRUE)
  } else {
    results <- rbindlist(mclapply(chunks, worker, mc.cores = n_cores),
                         use.names = TRUE, fill = TRUE)
  }
  if (nrow(results) == 0) return(results)

  results[, p_interact_fdr_global := p.adjust(p_interact, method = "BH")]
  results[, p_interact_fdr_locus := p.adjust(p_interact, method = "BH"),
          by = locus_id]
  results[, rank_global := frank(p_interact, ties.method = "min")]
  results[, rank_locus := frank(p_interact, ties.method = "min"), by = locus_id]
  setorder(results, p_interact)
  results
}

cat("Building candidate set...\n")
candidates <- build_candidates()
fwrite(candidates, file.path(OUTDIR, "targeted_cis_trans_candidate_variants_unmapped.csv"))

cat("Reading PLINK BIM and mapping candidates...\n")
bim <- fread(
  paste0(PLINK_PREFIX, ".bim"),
  header = FALSE,
  col.names = c("chr", "snp", "cm", "pos", "a1", "a2")
)
mapped <- map_candidates_to_bim(candidates, bim)
fwrite(mapped, file.path(OUTDIR, "targeted_cis_trans_candidate_variants_mapped.csv"))

cat(sprintf("Unmapped candidates: %d; mapped unique gene-variant candidates: %d\n",
            nrow(candidates), nrow(mapped)))

cat("Reading PLINK genotypes for mapped candidates...\n")
geno <- read_plink_dosage(PLINK_PREFIX, mapped$bed_index)

cat("Loading expression matrix...\n")
fpkm <- fread(FPKM_FILE)
fpkm[, gene_id_clean := sub("\\..*$", "", gene_id)]

tumor_expr <- make_expression_wide(fpkm, unique(mapped$gene_id_clean), tf_panel, "T")
normal_expr <- make_expression_wide(fpkm, unique(mapped$gene_id_clean), tf_panel, "N")

tumor_data <- prepare_tissue_data(tumor_expr, geno, mapped, "Tumor")
normal_data <- prepare_tissue_data(normal_expr, geno, mapped, "Normal")

cat("Running tumor screen...\n")
tumor_results <- run_tissue_screen(tumor_data, mapped, tf_panel)
cat("Running normal screen...\n")
normal_results <- run_tissue_screen(normal_data, mapped, tf_panel)

all_results <- rbindlist(list(tumor_results, normal_results), use.names = TRUE, fill = TRUE)
setorder(all_results, tissue, p_interact)

all_file <- file.path(OUTDIR, "targeted_cis_trans_interaction_screen_all.csv")
fwrite(all_results, all_file)

top_results <- all_results[, .SD[1:min(.N, 2000)], by = tissue]
fwrite(top_results, file.path(OUTDIR, "targeted_cis_trans_interaction_screen_top2000_by_tissue.csv"))

foxm1_results <- all_results[TF == "FOXM1"]
setorder(foxm1_results, tissue, p_interact)
fwrite(foxm1_results, file.path(OUTDIR, "targeted_cis_trans_interaction_screen_FOXM1.csv"))

kif18b_results <- all_results[gene_symbol == "KIF18B" |
                                gene_id == "ENSG00000186185" |
                                model_variant_id == "rs9914422" |
                                bed_snp == "rs9914422"]
setorder(kif18b_results, tissue, p_interact)
fwrite(kif18b_results, file.path(OUTDIR, "targeted_cis_trans_interaction_screen_KIF18B_rs9914422.csv"))

locus_best <- all_results[, .SD[1], by = .(tissue, locus_id)]
setorder(locus_best, tissue, p_interact)
fwrite(locus_best, file.path(OUTDIR, "targeted_cis_trans_interaction_screen_best_regulator_per_locus.csv"))

if (nrow(tumor_results) > 0) {
  plot_dt <- tumor_results[1:min(30, nrow(tumor_results))]
  plot_dt[, label := paste(gene_symbol, model_variant_id, TF, sep = " | ")]
  plot_dt[, highlight := ifelse(TF == "FOXM1" & model_variant_id == "rs9914422",
                                "FOXM1-rs9914422-KIF18B", "Other")]
  p <- ggplot(plot_dt, aes(x = reorder(label, -log10(p_interact)),
                           y = -log10(p_interact), fill = highlight)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c("FOXM1-rs9914422-KIF18B" = "#D62828",
                                 "Other" = "#3B6EA8")) +
    labs(title = "Targeted cis-trans interaction screen",
         subtitle = "Top 30 tumor variant-gene-regulator interactions",
         x = NULL, y = "-log10 interaction P value", fill = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "top")
  ggsave(file.path(OUTDIR, "Fig_targeted_cis_trans_screen_top30_tumor.pdf"),
         p, width = 13, height = 6)
}

summary_file <- file.path(OUTDIR, "targeted_cis_trans_interaction_screen_summary.txt")
sink(summary_file)
cat("Targeted cis-trans interaction screen for recurrent PCa-ASE variants\n")
cat("==================================================================\n\n")
cat(sprintf("Date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("Candidate definition:\n")
cat(sprintf("- recurrent significant ASE SNPs: n_sig_samples >= %d\n", MIN_SIG_SAMPLES))
cat("- existing LD coding proxies: best_proxy_r2 >= 0.8\n")
cat("- forced validated KIF18B proxy: rs9914422\n\n")
cat("Model:\n")
cat("log2(target_gene FPKM + 0.01) ~ dosage + regulator_z + dosage:regulator_z\n\n")
cat("Filters:\n")
cat(sprintf("- minimum complete samples: %d\n", MIN_N))
cat(sprintf("- minimum minor allele count: %d\n", MIN_MAC))
cat(sprintf("- minimum median target FPKM: %.2f\n", MIN_TARGET_MEDIAN))
cat(sprintf("- minimum median regulator FPKM: %.2f\n\n", MIN_REGULATOR_MEDIAN))

cat("Candidate counts:\n")
print(candidates[, .N, by = candidate_source])
cat(sprintf("\nMapped unique gene-variant candidates: %d\n", nrow(mapped)))
cat("Mapped candidates by source:\n")
print(mapped[, .N, by = candidate_source])

cat("\nTest counts:\n")
print(all_results[, .(n_tests = .N,
                      n_loci = uniqueN(locus_id),
                      n_regulators = uniqueN(TF)),
                  by = tissue])

cat("\nTop 20 tumor interactions:\n")
if (nrow(tumor_results) > 0) {
  print(tumor_results[1:min(20, .N),
                      .(rank_global, gene_symbol, model_variant_id, candidate_source,
                        TF, n, mac, beta_interact, p_interact,
                        p_interact_fdr_global, p_interact_fdr_locus,
                        rank_locus, delta_R2)],
        digits = 4)
}

cat("\nFOXM1-rs9914422-KIF18B rows:\n")
fox_anchor <- all_results[
  TF == "FOXM1" &
    (model_variant_id == "rs9914422" | bed_snp == "rs9914422") &
    (gene_symbol == "KIF18B" | gene_id == "ENSG00000186185")
]
if (nrow(fox_anchor) > 0) {
  print(fox_anchor[, .(tissue, rank_global, rank_locus, gene_symbol,
                       model_variant_id, bed_snp, dosage_allele, n, mac,
                       beta_interact, p_interact, p_interact_fdr_global,
                       p_interact_fdr_locus, R2_additive, R2_interaction,
                       delta_R2)],
        digits = 4)
} else {
  cat("Not tested after filters.\n")
}

cat("\nBest KIF18B/rs9914422 regulators:\n")
if (nrow(kif18b_results) > 0) {
  print(kif18b_results[1:min(20, .N),
                       .(tissue, rank_global, rank_locus, gene_symbol,
                         model_variant_id, TF, beta_interact, p_interact,
                         p_interact_fdr_global, p_interact_fdr_locus,
                         delta_R2)],
        digits = 4)
}
sink()

cat("\nScreen complete.\n")
cat(sprintf("Results saved to: %s\n", OUTDIR))
cat(sprintf("Summary: %s\n", summary_file))
