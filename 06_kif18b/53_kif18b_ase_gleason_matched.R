#!/usr/bin/env Rscript
# 53: is KIF18B tumour ASE associated with Gleason score beyond what any gene with similar expression behaviour shows?
# Motivation: in the 04 audit all 9 tumours with significant KIF18B ASE were Gleason 8-9 (logistic P=0.037; +log depth P=0.052),
# but KIF18B expression rises with Gleason (rho 0.53), so higher-grade tumours have more allelic reads and more power to call ASE.
# Design fixed 2026-10-06 before running (the KIF18B numbers above had been seen; the matched null had not).
# Data: Sen per-gene results of the 118 ASE-cohort tumours (T{ID}_ase_results.csv; evaluable = fdr not NA; ASE = fdr < 0.10, as 04),
#   tumour allelic depth (ase_gene_sample_agg.csv.gz), tumour FPKM, CPGEA clinical GS and BCR (YES/NO; "Dead"/NA excluded).
# Per gene (>= 50 evaluable tumours with GS, >= 4 with ASE):
#   PRIMARY  z_gs_dep = signed LRT z of GS in glm(ASE ~ GS + log10 depth, binomial) vs glm(ASE ~ log10 depth)
#   secondary z_gs (no depth term); Spearman |a| vs GS (continuous magnitude); z_bcr_dep = as primary with BCR instead of GS
# Expression-matched null: the NMATCH genes nearest to KIF18B (Euclidean on standardised features) in
#   median tumour log2 FPKM, Spearman(tumour FPKM, GS), n evaluable tumours, n tumours with ASE  (event count sets the power).
# Empirical two-sided P = (1 + #matched genes with |stat| >= |KIF18B|) / (NMATCH + 1); percentile among all eligible genes also reported.
suppressPackageStartupMessages({library(data.table); library(parallel)}); setDTthreads(4); set.seed(53)
NMATCH <- as.integer(Sys.getenv("NMATCH", "200")); NC <- as.integer(Sys.getenv("NC", "8"))
P <- "/path/to/data_root/ASE"; O <- file.path(P, "code/v8_round2/kif18b_audit/ase_gleason_matched"); dir.create(O, showWarnings = FALSE)
GID <- "ENSG00000186185"
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

fs <- list.files(P, pattern = "^T\\d+_ase_results\\.csv$", full.names = TRUE)
A <- rbindlist(lapply(fs, function(f) fread(f, select = c("gene_id", "a", "fdr"))[, pid := gsub("\\D", "", basename(f))]))
A <- A[!is.na(fdr)][, `:=`(ase = fdr < 0.10, gene_id = sub("\\..*", "", gene_id))]
ag <- fread(file.path(P, "code/ase_gene_sample_agg.csv.gz"), select = c("gene_id", "sample_id", "dep"))
ag <- ag[grepl("^T", sample_id)][, .(gene_id = sub("\\..*", "", gene_id), pid = gsub("\\D", "", sample_id), dep)]
A <- merge(A, ag, by = c("gene_id", "pid"))[dep > 0]
cl <- fread("/path/to/data_root/CPGEA/CPGEA_web/CPGEA_Clinical_Info.tsv")
cl <- cl[, .(pid = as.character(`Sequencing ID`), GS = suppressWarnings(as.numeric(GS)), BCR = fifelse(BCR == "YES", 1L, fifelse(BCR == "NO", 0L, NA_integer_)))]
A <- merge(A, cl, by = "pid")[is.finite(GS)]
msg("tumours %d, tumour-gene observations %s", uniqueN(A$pid), format(nrow(A), big.mark = ","))

fp <- fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls")
tc <- grep("^T\\d+_WTS$", names(fp), value = TRUE)
fm <- log2(as.matrix(fp[, ..tc]) + 1); rownames(fm) <- sub("\\..*", "", fp$gene_id)
gsT <- cl[match(gsub("\\D", "", tc), pid), GS]; okc <- is.finite(gsT)
feat <- data.table(gene_id = rownames(fm), mexp = apply(fm, 1, median),
                   rho_exp_gs = as.numeric(cor(t(fm[, okc]), gsT[okc], method = "spearman")))

lrz <- function(y, x, d) { f1 <- glm(y ~ x + log10(d), family = binomial); f0 <- glm(y ~ log10(d), family = binomial)
  lr <- max(f0$deviance - f1$deviance, 0); sign(coef(f1)[2]) * sqrt(lr) }
lrz0 <- function(y, x) { f1 <- glm(y ~ x, family = binomial); lr <- max(f1$null.deviance - f1$deviance, 0); sign(coef(f1)[2]) * sqrt(lr) }
cnt <- A[, .(n = .N, k = sum(ase)), by = gene_id]
elig <- cnt[n >= 50 & k >= 4]$gene_id
msg("eligible genes (>= 50 evaluable tumours with GS, >= 4 with ASE): %s ; KIF18B eligible: %s", format(length(elig), big.mark = ","), GID %in% elig)
setkey(A, gene_id)
one <- function(g) { d <- A[.(g)]
  zb <- NA_real_; db <- d[!is.na(BCR)]; if (nrow(db) >= 40 && sum(db$ase) >= 3) zb <- suppressWarnings(lrz(db$ase, db$BCR, db$dep))
  data.table(gene_id = g, n = nrow(d), k = sum(d$ase), med_dep = median(d$dep),
             z_gs_dep = suppressWarnings(lrz(d$ase, d$GS, d$dep)), z_gs = suppressWarnings(lrz0(d$ase, d$GS)),
             rho_absa_gs = suppressWarnings(cor(abs(d$a), d$GS, method = "spearman")), z_bcr_dep = zb) }
S <- rbindlist(mclapply(elig, function(g) tryCatch(one(g), error = function(e) NULL), mc.cores = NC))
S <- merge(S, feat, by = "gene_id")
fwrite(S, file.path(O, "per_gene_ase_gleason.csv"))

K <- S[gene_id == GID]; if (!nrow(K)) stop("KIF18B not among eligible genes")
F <- c("mexp", "rho_exp_gs", "n", "k"); Z <- scale(as.matrix(S[, ..F]))
dist <- sqrt(colSums((t(Z) - Z[S$gene_id == GID, ])^2))
M <- S[order(dist)][gene_id != GID][seq_len(min(NMATCH, .N - 1))]
msg("\nKIF18B: n %d, ASE %d, median depth %.0f, median log2 FPKM %.2f, rho(FPKM, GS) %+.2f", K$n, K$k, K$med_dep, K$mexp, K$rho_exp_gs)
msg("matched genes (n=%d): median log2 FPKM %.2f, rho(FPKM, GS) %+.2f, n %.0f, ASE %.0f  (all eligible: %.2f, %+.2f, %.0f, %.0f)",
    nrow(M), median(M$mexp), median(M$rho_exp_gs), median(M$n), median(M$k), median(S$mexp), median(S$rho_exp_gs), median(S$n), median(S$k))
R <- rbindlist(lapply(c("z_gs_dep", "z_gs", "rho_absa_gs", "z_bcr_dep"), function(v) { kv <- K[[v]]; mv <- M[[v]][is.finite(M[[v]])]; av <- S[[v]][is.finite(S[[v]])]
  data.table(statistic = v, primary = v == "z_gs_dep", KIF18B = kv, matched_median = median(mv), matched_n = length(mv),
             P_matched = (1 + sum(abs(mv) >= abs(kv))) / (length(mv) + 1), pct_matched = 100 * mean(mv < kv),
             all_median = median(av), pct_all = 100 * mean(av < kv)) }))
fwrite(R, file.path(O, "kif18b_vs_matched.csv")); fwrite(M, file.path(O, "matched_genes.csv"))
for (i in seq_len(nrow(R))) with(R[i], msg("%-12s KIF18B %+.3f | matched median %+.3f (n=%d), two-sided P=%.3f, percentile %.0f | all-gene median %+.3f, percentile %.0f%s",
  statistic, KIF18B, matched_median, matched_n, P_matched, pct_matched, all_median, pct_all, if (primary) "   <- PRIMARY" else ""))
msg("DONE")
