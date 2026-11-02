#!/usr/bin/env Rscript
# Cell-type interaction eQTL (Kim-Hellmuth 2020 Science, GTEx ieQTL design), pre-specified before looking at results (2026-10-06).
# Leads: GTEx v8 Prostate eGenes (qval < 0.05, SNV lead), i.e. independent of CPGEA. Per tissue (tumour, normal) and per gene:
#   INT(log2 FPKM + 1) ~ G + 5 expression PCs + f + G:f        f = BayesPrism v2d cell fraction (GSE185344 v2 reference)
# PRIMARY (one per tissue): f = epithelial fraction; statistic lambda = median(t_GxF^2) / 0.455 over lead eQTLs;
#   null = whole-sample permutation of the fraction vector (main effect and interaction together), 100 permutations, one-sided P.
# Secondary / descriptive: number of genes with BH-FDR < 0.05; mean interaction t signed by the eQTL main-effect direction
#   (positive = eQTL stronger in samples with more of that cell type).
# Sensitivity: f = fibroblast, f = smooth muscle; tumour model + purity; CPGEA normal-map leads.
suppressPackageStartupMessages({library(data.table); library(parallel)}); setDTthreads(1); set.seed(48)
V <- "/path/to/data_root/ASE/code/v8_round2"; O <- file.path(V, "celltype_ieqtl"); dir.create(O, showWarnings = FALSE)
QD <- "/path/to/data_root/CPGEA/QTLs/CPGEA_eQTL_and_sQTL/"; PP <- "/path/to/data_root/CPGEA/CPGEA_genotype/combined.134.qc"
NREP <- as.integer(Sys.getenv("NREP", "100")); NC <- as.integer(Sys.getenv("NC", "16"))
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

bim <- fread(paste0(PP, ".bim"), col.names = c("chr","snp","cm","pos","a1","a2"))[, idx := .I]
fam <- fread(paste0(PP, ".fam"), col.names = c("fid","iid","f","m","sex","pheno"))
sids <- as.character(sub("^N", "", fam$fid)); N <- nrow(fam); bps <- ceiling(N / 4)
CON <- file(paste0(PP, ".bed"), "rb"); readBin(CON, "raw", n = 3)
LUT <- matrix(NA_integer_, 256L, 4L)
for (bv in 0:255) for (j in 0:3) { cd <- bitwAnd(bitwShiftR(bv, 2L*j), 3L); LUT[bv+1L, j+1L] <- if (cd == 0L) 2L else if (cd == 1L) NA_integer_ else if (cd == 2L) 1L else 0L }
readsnp <- function(idx) { seek(CON, 3 + (idx - 1) * bps); rb <- readBin(CON, "raw", n = bps); as.vector(t(LUT[as.integer(rb) + 1L, , drop = FALSE]))[seq_len(N)] }
b2 <- bim[, .(idx = idx[1], a1 = a1[1], a2 = a2[1], nd = .N), by = .(chr = as.character(chr), pos)][nd == 1]

fpkm <- fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls"); fpkm[, eg := sub("\\..*", "", gene_id)]
pur <- fread("/path/to/data_root/variants/batch_43/facets/Table_S1_PRAD.txt")[, .(sid = as.character(SAMPLE_ID), purity = as.numeric(purity))]
int <- function(y) qnorm((rank(y) - 0.5) / length(y))

gt <- fread(file.path(V, "delta_eqtl/Prostate.v8.egenes.txt.gz"))[qval < 0.05]
gt[, c("vchr", "vpos", "vref", "valt") := tstrsplit(variant_id, "_", keep = 1:4)]
LG <- gt[nchar(vref) == 1 & nchar(valt) == 1, .(eg = sub("\\..*", "", gene_id), chr = sub("^chr", "", vchr), pos = as.integer(vpos), ref = vref, alt = valt)]
LG <- merge(LG, b2, by = c("chr", "pos"))[(a1 == ref & a2 == alt) | (a1 == alt & a2 == ref)]
leadN <- function() { q <- fread(paste0(QD, "CPGEA_cis_normal_eQTLs_anno.txt"))[FDR < 0.05]
  q[, eg := sub("\\..*", "", sub("\\|.*", "", gene))]; q <- q[order(p.value)][!duplicated(eg)][, .(eg, chr = sub("^chr", "", Chr), pos = Pos)]
  merge(q, b2, by = c("chr", "pos")) }
LC <- leadN()
msg("GTEx prostate SNV leads matched to CPGEA genotypes: %d ; CPGEA normal-map leads: %d", nrow(LG), nrow(LC))

run <- function(tis, ctype, L, lead_lab, add_purity = FALSE) {
  cc <- grep(sprintf("^%s\\d+_WTS$", tis), names(fpkm), value = TRUE); es <- gsub("\\D", "", cc); keep <- es %in% sids; cc <- cc[keep]; es <- es[keep]
  dec <- fread(file.path(V, sprintf("decon_v2d/CPGEA_%s_cell_fraction.tsv", if (tis == "T") "tumor" else "normal")))
  dec <- dec[, .(sid = gsub("\\D", "", get(names(dec)[1])), f = get(ctype))]
  cv <- merge(data.table(sid = es), dec, by = "sid", all.x = TRUE, sort = FALSE)
  if (add_purity) cv <- merge(cv, pur, by = "sid", all.x = TRUE, sort = FALSE)
  cv <- cv[match(es, sid)]; ok <- complete.cases(cv); cc <- cc[ok]; es <- es[ok]; cv <- cv[ok]; n <- length(es)
  E <- log2(as.matrix(fpkm[, ..cc]) + 1); rownames(E) <- fpkm$eg
  ex <- E[rowMeans(E > 0.5) > 0.5, , drop = FALSE]; pcs <- prcomp(t(ex - rowMeans(ex)), rank. = 5)$x
  base <- cbind(1, pcs, if (add_purity) cbind(purity = cv$purity) else NULL)
  L <- L[eg %in% rownames(ex)][!duplicated(eg)]
  G <- sapply(L$idx, readsnp); rownames(G) <- sids; G <- G[es, , drop = FALSE]
  G <- apply(G, 2, function(g) { g[is.na(g)] <- mean(g, na.rm = TRUE); g }); ok <- apply(G, 2, function(g) sum(round(g) != 0) >= 5 & sum(round(g) != 2) >= 5)
  G <- G[, ok, drop = FALSE]; L <- L[ok]; Y <- apply(E[L$eg, , drop = FALSE], 1, int)
  fz <- as.numeric(scale(cv$f))
  fit <- function(f) { tt <- sg <- numeric(ncol(G))
    for (j in seq_len(ncol(G))) { X <- cbind(base, f = f, G = G[, j], GxF = G[, j] * f); r <- lm.fit(X, Y[, j])
      k <- ncol(X); if (r$rank < k) { tt[j] <- NA; next }
      df <- n - k; s2 <- sum(r$residuals^2) / df; XtXi <- chol2inv(qr.R(r$qr))
      tt[j] <- r$coefficients[k] / sqrt(s2 * XtXi[k, k]); sg[j] <- sign(r$coefficients[k - 1]) }
    list(t = tt, sg = sg, df = n - ncol(base) - 3) }
  ob <- fit(fz); lam <- median(ob$t^2, na.rm = TRUE) / 0.455; pv <- 2 * pt(-abs(ob$t), ob$df); fdr <- p.adjust(pv, "BH")
  nl <- unlist(mclapply(seq_len(NREP), function(b) { set.seed(1000 + b); median(fit(fz[sample.int(n)])$t^2, na.rm = TRUE) / 0.455 }, mc.cores = NC))
  st <- (ob$t * ob$sg)[is.finite(ob$t)]; gb <- quantile(replicate(2000, mean(sample(st, replace = TRUE))), c(.025, .975))
  r <- data.table(tissue = tis, cell = ctype, leads = lead_lab, purity = add_purity, n = n, genes = ncol(G), frac_median = median(cv$f), frac_iqr = IQR(cv$f),
                  lambda = lam, null_lambda = mean(nl), null_lambda_max = max(nl), P_lambda = (1 + sum(nl >= lam)) / (NREP + 1),
                  fdr05 = sum(fdr < 0.05), signed_mean_t = mean(st), signed_lo = gb[1], signed_hi = gb[2])
  msg("[%s] %-16s %-11s%s n=%d genes=%d | %s median %.2f (IQR %.2f) | lambda %.3f (null %.3f, max %.3f, P=%.3f) | FDR<0.05 %d | signed t %+.3f [%+.3f, %+.3f]",
      tis, ctype, lead_lab, if (add_purity) " +purity" else "", n, ncol(G), ctype, median(cv$f), IQR(cv$f), lam, mean(nl), max(nl), r$P_lambda, r$fdr05, mean(st), gb[1], gb[2])
  fwrite(data.table(L[, .(eg, chr, pos)], t = ob$t, signed_t = ob$t * ob$sg, p = pv, fdr = fdr)[order(p)],
         file.path(O, sprintf("ieqtl_%s_%s_%s%s.csv.gz", tis, ctype, lead_lab, if (add_purity) "_pur" else "")))
  r }

R <- list()
msg("\n== PRIMARY: epithelial fraction, GTEx leads ==")
for (tis in c("T", "N")) R[[length(R) + 1]] <- run(tis, "Epithelial_cells", LG, "GTEx")
if (Sys.getenv("ONLY_PRIMARY") == "1") { R <- rbindlist(R); print(R); fwrite(R, file.path(O, sprintf("celltype_ieqtl_primary_%dperm.csv", NREP))); msg("DONE"); quit(save = "no") }
msg("\n== SENSITIVITY ==")
for (tis in c("T", "N")) for (ct in c("Fibroblast", "SMC")) R[[length(R) + 1]] <- run(tis, ct, LG, "GTEx")
R[[length(R) + 1]] <- run("T", "Epithelial_cells", LG, "GTEx", add_purity = TRUE)
for (tis in c("T", "N")) R[[length(R) + 1]] <- run(tis, "Epithelial_cells", LC, "CPGEA_N")
R <- rbindlist(R); print(R); fwrite(R, file.path(O, "celltype_ieqtl_summary.csv")); msg("DONE")
