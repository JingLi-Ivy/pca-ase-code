#!/usr/bin/env Rscript
# 49: does the original cis x trans statistic (difference of Spearman correlations, het vs hom) carry a reliability bias?
# Part 1  On the gene/variant selection code of 13/26b (seed 42), with normal-tissue depth from the normal file as in 26e:
#         per gene, heterozygous vs homozygous group: median reporter SNPs (nsnp), median aggregate depth (dep),
#         approximate reliability of |a| = 1 - mean(binomial noise variance 0.25/dep) / var(|a|).
# Part 2  Simulation on the real structure (real het/hom split, real depth, real TF z-scores; 1,000 genes x 20 TFs):
#         true imbalance a = 0.10 + g * TF + N(0, 0.05), the SAME slope g in het and hom; reads ~ Binomial(dep, 0.5 + a);
#         |a_obs| residualised on log10 depth as in 26b (r0). Statistics: original Fisher-z of Spearman difference (13/noprior)
#         and slope z (26b). True correlation r in {-0.10, 0, +0.05, +0.10, +0.20}. Binomial noise ignores overdispersion,
#         so real noise is larger and the simulated bias is a lower bound.
suppressPackageStartupMessages({library(data.table); library(parallel)}); setDTthreads(1)
NC <- as.integer(Sys.getenv("NC", "16")); NSIM <- as.integer(Sys.getenv("NSIM", "50"))
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
OUT <- "/path/to/data_root/ASE/code/v8_round2/reliability_bias"; dir.create(OUT, showWarnings = FALSE)
PP <- "/path/to/data_root/CPGEA/CPGEA_genotype/combined.134.qc"; WIN <- 10000; set.seed(42)

bim <- fread(paste0(PP, ".bim"), col.names = c("chr","snp","cm","pos","a1","a2"))[, idx := .I][, chr := paste0("chr", chr)]
fam <- fread(paste0(PP, ".fam"), col.names = c("fid","iid","f","m","sex","pheno"))
sids <- as.character(sub("^N", "", fam$fid)); N <- nrow(fam); bps <- ceiling(N / 4)
CON <- file(paste0(PP, ".bed"), "rb"); readBin(CON, "raw", n = 3)
LUT <- matrix(NA_integer_, 256L, 4L)
for (bv in 0:255) for (j in 0:3) { cd <- bitwAnd(bitwShiftR(bv, 2L*j), 3L); LUT[bv+1L, j+1L] <- if (cd == 0L) 2L else if (cd == 1L) NA_integer_ else if (cd == 2L) 1L else 0L }
readsnp <- function(idx) { seek(CON, 3 + (idx - 1) * bps); rb <- readBin(CON, "raw", n = bps); as.vector(t(LUT[as.integer(rb) + 1L, , drop = FALSE]))[seq_len(N)] }
ann <- unique(fread("/path/to/data_root/ASE/all_T_samples_SNPs_with_ASE_annotation.csv", select = c("seqtnames", "start", "gene_id")))
ann[, chr := paste0("chr", as.character(seqtnames))]
gpos <- ann[, .(gchr = chr[1], gmin = min(start), gmax = max(start)), by = gene_id]
fpkm <- fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls")
tcc <- grep("^T\\d+_WTS$", names(fpkm), value = TRUE); ncc <- grep("^N\\d+_WTS$", names(fpkm), value = TRUE)
esT <- as.character(gsub("\\D", "", tcc)); esN <- as.character(gsub("\\D", "", ncc)); PAIR <- intersect(esT, esN)
ag <- fread("/path/to/data_root/ASE/code/ase_gene_sample_agg.csv.gz")[, sid := as.character(gsub("\\D", "", sample_id))]
agN <- fread("/path/to/data_root/ASE/code/cistrans_full/ase_gene_sample_agg_N.csv.gz")[, sid := as.character(sid)]
loadase <- function(pre, es) {
  fs <- list.files("/path/to/data_root/ASE", pattern = sprintf("^%s\\d+_ase_results\\.csv$", pre), full.names = TRUE)
  a <- rbindlist(lapply(fs, function(f) { d <- fread(f)[, .(gene_id, a)]; d[, sid := as.character(gsub("\\D", "", basename(f)))]; d }))
  d <- merge(a, (if (pre == "N") agN else ag)[, .(gene_id, sid, dep, nsnp)], by = c("gene_id", "sid"))[dep >= 20 & sid %in% es]
  d[, `:=`(ai = abs(a), ld = log10(dep))] }
aT <- loadase("T", esT)[sid %in% PAIR]; aN <- loadase("N", esN)[sid %in% PAIR]
gok <- intersect(aT[, .N, by = gene_id][N >= 80]$gene_id, aN[, .N, by = gene_id][N >= 80]$gene_id)
aT <- aT[gene_id %in% gok]; aN <- aN[gene_id %in% gok]
setorder(aT, gene_id); setorder(aN, gene_id)
aT[, r0 := residuals(lm(ai ~ ld)), by = gene_id]; aN[, r0 := residuals(lm(ai ~ ld)), by = gene_id]
setkey(aT, gene_id); setkey(aN, gene_id)
both <- intersect(esT, esN); setkey(bim, chr); sel <- list()
for (g in gok) {
  gp <- gpos[gene_id == g]; if (!nrow(gp)) next
  cd <- bim[.(gp$gchr)][pos >= gp$gmin - WIN & pos <= gp$gmax + WIN]; if (!nrow(cd)) next
  if (nrow(cd) > 40) cd <- cd[sample(.N, 40)]
  best <- NULL; bv <- -1L
  for (k in seq_len(nrow(cd))) { H <- setNames(readsnp(cd$idx[k]), sids)[both]
    v <- min(sum(H == 1, na.rm = TRUE), sum(H %in% c(0, 2), na.rm = TRUE)); if (v > bv) { bv <- v; best <- cd[k] } }
  if (is.null(best) || bv < 25) next
  sel[[length(sel) + 1]] <- list(gene_id = g, idx = best$idx) }
S <- rbindlist(sel); GT <- sapply(S$idx, readsnp); rownames(GT) <- sids
msg("genes in both tissues %s ; selected genes %s (26b: 10,196 with tumour depth for normal; this script uses normal depth, as 26e)", format(length(gok), big.mark = ","), format(nrow(S), big.mark = ","))

grp <- function(A, tis) rbindlist(lapply(seq_len(nrow(S)), function(i) { z <- A[.(S$gene_id[i])]; hh <- GT[z$sid, i]; z <- z[!is.na(hh)][, het := hh[!is.na(hh)] == 1L]
  if (z[het == TRUE, .N] < 25 || z[het == FALSE, .N] < 25) return(NULL)
  z[, .(n = .N, nsnp = as.double(median(nsnp)), dep = as.double(median(dep)), rel = 1 - mean(0.25 / dep) / var(ai)), by = het][, `:=`(gene = S$gene_id[i], tissue = tis)] }))
P1 <- rbind(grp(aT, "tumour"), grp(aN, "normal")); fwrite(P1, file.path(OUT, "het_hom_snp_depth.csv.gz"))
W <- dcast(P1, tissue + gene ~ het, value.var = c("n", "nsnp", "dep", "rel"))
msg("\n== Part 1: heterozygous vs homozygous group at the selected variant (per-gene medians; summarised across genes) ==")
for (tis in c("tumour", "normal")) { w <- W[tissue == tis]
  msg("[%s] genes %s | patients per gene het %.0f / hom %.0f", tis, format(nrow(w), big.mark = ","), median(w$n_TRUE), median(w$n_FALSE))
  msg("   reporter SNPs per patient-gene: median het %.1f vs hom %.1f ; mean het %.2f vs hom %.2f ; genes het > hom %.1f%% (het < hom %.1f%%) ; Wilcoxon P=%.1g",
      median(w$nsnp_TRUE), median(w$nsnp_FALSE), mean(w$nsnp_TRUE), mean(w$nsnp_FALSE), 100 * mean(w$nsnp_TRUE > w$nsnp_FALSE), 100 * mean(w$nsnp_TRUE < w$nsnp_FALSE), wilcox.test(w$nsnp_TRUE, w$nsnp_FALSE, paired = TRUE)$p.value)
  msg("   aggregate depth per patient-gene: median het %.0f vs hom %.0f ; genes het > hom %.1f%% ; median ratio het/hom %.2f ; Wilcoxon P=%.1g",
      median(w$dep_TRUE), median(w$dep_FALSE), 100 * mean(w$dep_TRUE > w$dep_FALSE), median(w$dep_TRUE / w$dep_FALSE), wilcox.test(w$dep_TRUE, w$dep_FALSE, paired = TRUE)$p.value)
  msg("   approx. reliability of |a| (binomial noise only): median het %.3f vs hom %.3f ; genes het > hom %.1f%% ; implied attenuation ratio sqrt(rel_hom/rel_het) median %.3f",
      median(w$rel_TRUE), median(w$rel_FALSE), 100 * mean(w$rel_TRUE > w$rel_FALSE), median(sqrt(pmax(w$rel_FALSE, 0.01) / pmax(w$rel_TRUE, 0.01)))) }

fpk <- function(cc, es) { tf <- fread("/path/to/data_root/ASE/code/adastra_tf_map.csv")$sym; m <- fpkm[gene_name %in% tf, c("gene_name", cc), with = FALSE]
  m <- m[apply(as.matrix(m[, -1]), 1, median) >= 1][!duplicated(gene_name)]; z <- t(apply(log2(as.matrix(m[, -1]) + 0.01), 1, scale)); dimnames(z) <- list(m$gene_name, es); z }
TZ <- list(tumour = fpk(tcc, esT), normal = fpk(ncc, esN)); A <- list(tumour = aT, normal = aN)
int <- function(y) qnorm((rank(y) - 0.5) / length(y))
fz <- function(y1, x1, y0, x0) { c1 <- cor(y1, x1, method = "spearman"); c0 <- cor(y0, x0, method = "spearman")
  (atanh(c1) - atanh(c0)) / sqrt(1 / (length(y1) - 3) + 1 / (length(y0) - 3)) }
slz <- function(y1, x1, y0, x0) { f <- function(y, x) { b <- cov(y, x) / var(x); r <- y - mean(y) - b * (x - mean(x)); c(b, sum(r^2) / (length(y) - 2) / sum((x - mean(x))^2)) }
  a <- f(y1, x1); b <- f(y0, x0); (a[1] - b[1]) / sqrt(a[2] + b[2]) }
set.seed(49); SB <- 0.05
sim1 <- function(rr, tis, rep) { set.seed(rep * 1000 + round(rr * 100) + if (tis == "tumour") 0 else 500)
  Ax <- A[[tis]]; Tz <- TZ[[tis]]; ii <- sample(nrow(S), min(1000, nrow(S))); g_true <- rr * SB / sqrt(1 - rr^2); out <- list()
  for (i in ii) { z <- Ax[.(S$gene_id[i])]; hh <- GT[z$sid, i]; k <- !is.na(hh); z <- z[k]; het <- hh[k] == 1L
    if (sum(het) < 25 || sum(!het) < 25) next
    for (tf in sample(rownames(Tz), 20)) { x <- Tz[tf, z$sid]
      at <- pmin(pmax(0.10 + g_true * x + rnorm(length(x), 0, SB), -0.45), 0.45)
      ao <- abs(rbinom(length(x), z$dep, 0.5 + at) / z$dep - 0.5); y <- residuals(lm(ao ~ log10(z$dep)))
      yi <- int(y); out[[length(out) + 1]] <- c(fz = fz(y[het], x[het], y[!het], x[!het]), slz = slz(yi[het], x[het], yi[!het], x[!het])) } }
  m <- do.call(rbind, out); data.table(r = rr, tissue = tis, rep = rep, tests = nrow(m), fz = mean(m[, "fz"]), slz = mean(m[, "slz"])) }
grid <- CJ(r = c(-0.10, 0, 0.05, 0.10, 0.20), tissue = c("tumour", "normal"), rep = seq_len(NSIM))
SIM <- rbindlist(mclapply(seq_len(nrow(grid)), function(k) sim1(grid$r[k], grid$tissue[k], grid$rep[k]), mc.cores = NC))
fwrite(SIM, file.path(OUT, "simulation_reps.csv"))
msg("\n== Part 2: simulation, same true slope in het and hom (%d replicates x up to 1,000 genes x 20 TFs; true |a| SD %.2f) ==", NSIM, SB)
print(SIM[, .(tests = round(mean(tests)), fz_mean = mean(fz), fz_lo = quantile(fz, .025), fz_hi = quantile(fz, .975),
              slope_mean = mean(slz), slope_lo = quantile(slz, .025), slope_hi = quantile(slz, .975)), by = .(tissue, r)][order(tissue, r)])
msg("observed (26b, depth only; normal side used tumour depth): tumour correlation-difference +0.0204 / slope +0.0038 ; normal -0.012 / -0.0043. Compare with 26e.")
msg("DONE")
