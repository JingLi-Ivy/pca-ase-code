#!/usr/bin/env Rscript
# 50: expected H5 under NO cis x trans, given the reliability bias of the slope statistic on |a| (see 49).
# (1) Data: for every selected gene (13/26b selection code, seed 42; normal depth from the normal file as in 26e) and 30 random TFs,
#     the pooled (het + hom) Spearman correlation of the depth-residualised |a| with the TF -- the observable that drives the bias --
#     and the slope z of 26b/26e on the same tests.
# (2) Simulation on the real structure as in 49 (same true slope in het and hom), grid of true r and three true-|a| SDs (0.03/0.05/0.08);
#     records the same observable (pooled Spearman), the slope z, and |a| reliability.
# (3) Calibration: per tissue and SD, linear map slope-z ~ pooled Spearman from the simulation, evaluated at the observed pooled
#     Spearman -> expected bias; expected H5 = bias(tumour) - bias(normal). The SD whose simulated reliability matches the data
#     (49: tumour het/hom 0.80/0.69, normal 0.78/0.68) is the primary setting; the others are sensitivity.
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

fpk <- function(cc, es) { tf <- fread("/path/to/data_root/ASE/code/adastra_tf_map.csv")$sym; m <- fpkm[gene_name %in% tf, c("gene_name", cc), with = FALSE]
  m <- m[apply(as.matrix(m[, -1]), 1, median) >= 1][!duplicated(gene_name)]; z <- t(apply(log2(as.matrix(m[, -1]) + 0.01), 1, scale)); dimnames(z) <- list(m$gene_name, es); z }
TZ <- list(tumour = fpk(tcc, esT), normal = fpk(ncc, esN)); A <- list(tumour = aT, normal = aN)
int <- function(y) qnorm((rank(y) - 0.5) / length(y))
slz <- function(y1, x1, y0, x0) { f <- function(y, x) { b <- cov(y, x) / var(x); r <- y - mean(y) - b * (x - mean(x)); c(b, sum(r^2) / (length(y) - 2) / sum((x - mean(x))^2)) }
  a <- f(y1, x1); b <- f(y0, x0); (a[1] - b[1]) / sqrt(a[2] + b[2]) }

set.seed(50); OBS <- list()
for (tis in c("tumour", "normal")) { Ax <- A[[tis]]; Tz <- TZ[[tis]]; out <- list()
  for (i in seq_len(nrow(S))) { z <- Ax[.(S$gene_id[i])]; hh <- GT[z$sid, i]; k <- !is.na(hh); z <- z[k]; het <- hh[k] == 1L
    if (sum(het) < 25 || sum(!het) < 25) next
    y <- z$r0; yi <- int(y)
    for (tf in sample(rownames(Tz), 30)) { x <- Tz[tf, z$sid]
      out[[length(out) + 1]] <- c(i = i, rho = cor(y, x, method = "spearman"), slz = slz(yi[het], x[het], yi[!het], x[!het])) } }
  m <- as.data.table(do.call(rbind, out)); gm <- m[, .(rho = mean(rho), slz = mean(slz)), by = i]
  OBS[[tis]] <- data.table(tissue = tis, tests = nrow(m), genes = nrow(gm), rho = mean(gm$rho), slz = mean(gm$slz),
                           rho_lo = quantile(replicate(2000, mean(sample(gm$rho, replace = TRUE))), .025), rho_hi = quantile(replicate(2000, mean(sample(gm$rho, replace = TRUE))), .975)) }
OBS <- rbindlist(OBS); msg("\n(1) observed, 30 random TFs per gene (gene means):"); print(OBS)

sim1 <- function(rr, sb, tis, rep) { set.seed(rep * 10007 + round(rr * 1000) * 31 + round(sb * 1000) * 7 + if (tis == "tumour") 0 else 3)
  Ax <- A[[tis]]; Tz <- TZ[[tis]]; ii <- sample(nrow(S), min(1000, nrow(S))); g <- rr * sb / sqrt(1 - rr^2); out <- list()
  for (i in ii) { z <- Ax[.(S$gene_id[i])]; hh <- GT[z$sid, i]; k <- !is.na(hh); z <- z[k]; het <- hh[k] == 1L
    if (sum(het) < 25 || sum(!het) < 25) next
    for (tf in sample(rownames(Tz), 20)) { x <- Tz[tf, z$sid]
      at <- pmin(pmax(0.10 + g * x + rnorm(length(x), 0, sb), -0.45), 0.45)
      ao <- abs(rbinom(length(x), z$dep, 0.5 + at) / z$dep - 0.5); y <- residuals(lm(ao ~ log10(z$dep))); yi <- int(y)
      out[[length(out) + 1]] <- c(rho = cor(y, x, method = "spearman"), slz = slz(yi[het], x[het], yi[!het], x[!het]),
        relH = 1 - mean(0.25 / z$dep[het]) / var(ao[het]), relO = 1 - mean(0.25 / z$dep[!het]) / var(ao[!het])) } }
  m <- do.call(rbind, out); data.table(r = rr, sb = sb, tissue = tis, rep = rep, rho = mean(m[, "rho"]), slz = mean(m[, "slz"]),
                                        relH = median(m[, "relH"]), relO = median(m[, "relO"])) }
grid <- CJ(r = c(-0.10, -0.05, 0, 0.05, 0.10, 0.15), sb = c(0.03, 0.05, 0.08), tissue = c("tumour", "normal"), rep = 1:20)
SIM <- rbindlist(mclapply(seq_len(nrow(grid)), function(k) sim1(grid$r[k], grid$sb[k], grid$tissue[k], grid$rep[k]), mc.cores = NC))
fwrite(SIM, file.path(OUT, "h5_bias_simulation.csv"))
msg("\n(2) simulation summary (means over 20 replicates):")
print(SIM[, .(rho = mean(rho), slz = mean(slz), relH = mean(relH), relO = mean(relO)), by = .(sb, tissue, r)][order(sb, tissue, r)])

msg("\n(3) expected slope-z bias at the observed pooled Spearman (linear map from simulation; 95%% range over replicate-level fits)")
E <- rbindlist(lapply(c(0.03, 0.05, 0.08), function(SBv) rbindlist(lapply(c("tumour", "normal"), function(tis) {
  d <- SIM[sb == SBv & tissue == tis]; f <- lm(slz ~ rho, d); x0 <- OBS[tissue == tis, rho]
  bb <- replicate(2000, { q <- d[sample(.N, replace = TRUE)]; sum(coef(lm(slz ~ rho, q)) * c(1, x0)) })
  data.table(sb = SBv, tissue = tis, obs_rho = x0, slope_per_rho = coef(f)[2], bias = sum(coef(f) * c(1, x0)), bias_lo = quantile(bb, .025), bias_hi = quantile(bb, .975),
             sim_relH = mean(d$relH), sim_relO = mean(d$relO), obs_slz = OBS[tissue == tis, slz]) }))))
print(E)
X <- dcast(E, sb ~ tissue, value.var = c("bias", "obs_slz"))
X[, `:=`(expected_H5 = bias_tumour - bias_normal, observed_H5_here = obs_slz_tumour - obs_slz_normal)]
msg("\nexpected H5 under no cis x trans vs observed (same 30-TF subsample); full-TF observed H5: 26b +0.0104 (tumour depth for normal), 26e see its log")
print(X[, .(sb, bias_tumour, bias_normal, expected_H5, observed_H5_here)])
f26e <- "/path/to/data_root/ASE/code/v8_round2/cistrans_slope_26e/tn_slope_by_covariate_model.csv"
if (file.exists(f26e)) print(fread(f26e)[model == "ai ~ ld"])
msg("DONE")
