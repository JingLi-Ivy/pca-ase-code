#!/usr/bin/env Rscript
# 28f (2026-10-12): 28e restricted to PAIRED patients in both tissues (tumour and normal RNA, BayesPrism epithelial fraction, genotype; 118 patients, purity of 233 imputed by median).
# 28e (2026-10-12): 28d with GTEx v8 Prostate lead eQTLs instead of the CPGEA tissue-specific maps (independent lead set, as Fig4A/34c and S10/48).
# 28d (2026-10-08): 28c with 1,000 random non-TF draws (was 20; TFs were above all 20 draws) and the draw lambdas saved; output to cistrans_validate_28d.
# 28c: specificity controls for M1 (copied set-up of 28a). (1) random expressed non-TF genes in place of TFs (same number,
#   20 draws): if lambda is inflated equally, the M1 inflation reflects genotype x cell state, not TF dosage.
# (2) ASB at the lead eQTL: for gene x TF pairs where ADASTRA shows significant allele-specific binding of that TF at the lead
#   SNP, the interaction on the a1 dosage should be signed by the preferred allele (TF up -> preferred allele more expressed);
#   compare with the same TFs at other genes and with a label shuffle.
# Original M1 header:
# cis x trans, method M1: TF-interaction eQTL on total expression (GTEx ieQTL / Kim-Hellmuth 2020; Flynn 2022 validation).
# Per tissue, lead cis-eQTL per eGene (CPGEA eQTL map, FDR < 0.05, minimum P). For every eGene x expressed TF:
#   INT(log2 FPKM) ~ G + covariates + TF + G:TF        (M1a covariates: 5 expression PCs, epithelial fraction [, purity in tumour])
#                                                      (M1b: M1a + G:epithelial, the standard guard against cell-type ieQTLs)
# Statistic: interaction t signed by the direction of the eQTL main effect (positive = the TF amplifies the eQTL).
# Summaries: mean signed t, lambda = median(t^2)/0.455, fraction |t| > 3. Strict null: rows of the TF matrix permuted
# (one patient permutation for all genes and TFs; keeps TF co-expression, genotype, covariates and LD).
# Everything is computed by QR residualisation, vectorised over TFs.
suppressPackageStartupMessages({library(data.table); library(parallel)}); setDTthreads(1); set.seed(28)
O <- "/path/to/data_root/ASE/code/v8_round2/cistrans_validate_28f"; NDRAW <- as.integer(Sys.getenv("NDRAW", "1000")); dir.create(O, showWarnings = FALSE)
V <- "/path/to/data_root/ASE/code/v8_round2"; QD <- "/path/to/data_root/CPGEA/QTLs/CPGEA_eQTL_and_sQTL/"
PP <- "/path/to/data_root/CPGEA/CPGEA_genotype/combined.134.qc"; AD <- "/path/to/data_root/ASE/ADASTRA6/ADASTRA.v.6.1.Mabel/TF"
NREP <- as.integer(Sys.getenv("NREP", "100")); NC <- as.integer(Sys.getenv("NC", "16"))
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

bim <- fread(paste0(PP, ".bim"), col.names = c("chr","snp","cm","pos","a1","a2"))[, idx := .I]
fam <- fread(paste0(PP, ".fam"), col.names = c("fid","iid","f","m","sex","pheno"))
sids <- as.character(sub("^N", "", fam$fid)); N <- nrow(fam); bps <- ceiling(N / 4)
CON <- file(paste0(PP, ".bed"), "rb"); readBin(CON, "raw", n = 3)
LUT <- matrix(NA_integer_, 256L, 4L)
for (bv in 0:255) for (j in 0:3) { cd <- bitwAnd(bitwShiftR(bv, 2L*j), 3L); LUT[bv+1L, j+1L] <- if (cd == 0L) 2L else if (cd == 1L) NA_integer_ else if (cd == 2L) 1L else 0L }
readsnp <- function(idx) { seek(CON, 3 + (idx - 1) * bps); rb <- readBin(CON, "raw", n = bps); as.vector(t(LUT[as.integer(rb) + 1L, , drop = FALSE]))[seq_len(N)] }

fpkm <- fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls"); fpkm[, eg := sub("\\..*", "", gene_id)]
map <- fread("/path/to/data_root/ASE/code/adastra_tf_map.csv")[file.exists(file.path(AD, paste0(mn, "_HUMAN.tsv")))]
pur <- fread("/path/to/data_root/variants/batch_43/facets/Table_S1_PRAD.txt")[, .(sid = as.character(SAMPLE_ID), purity = as.numeric(purity))]
int <- function(y) qnorm((rank(y) - 0.5) / length(y))
GL <- fread(file.path(V, "delta_eqtl/Prostate.v8.egenes.txt.gz"))[qval < 0.05]
GL[, c("vchr", "vpos", "vref", "valt") := tstrsplit(variant_id, "_", keep = 1:4)]
GL <- GL[nchar(vref) == 1 & nchar(valt) == 1, .(eg = sub("\\..*", "", gene_id), chr = sub("^chr", "", vchr), pos = as.integer(vpos), ref = vref, alt = valt)]

dpid <- function(tis) { d <- fread(file.path(V, sprintf("decon_v2/CPGEA_%s_cell_fraction.tsv", if (tis == "T") "tumor" else "normal")))
  d[is.finite(Epithelial_cells), gsub("\\D", "", get(names(d)[1]))] }
PAIRIDS <- Reduce(intersect, list(gsub("\\D", "", grep("^T\\d+_WTS$", names(fpkm), value = TRUE)), gsub("\\D", "", grep("^N\\d+_WTS$", names(fpkm), value = TRUE)),
  dpid("T"), dpid("N"), sids,
  gsub("\\D", "", list.files("/path/to/data_root/ASE", pattern = "^T\\d+_ase_results\\.csv$"))))
stopifnot(length(PAIRIDS) == 118L)
pur[!is.finite(purity), purity := median(pur$purity, na.rm = TRUE)]
msg("purity imputed (cohort median %.3f) for: %s", median(pur$purity), paste(setdiff(PAIRIDS, fread("/path/to/data_root/variants/batch_43/facets/Table_S1_PRAD.txt")[is.finite(as.numeric(purity)), as.character(SAMPLE_ID)]), collapse = ","))
msg("paired patients used in both tissues: %d", length(PAIRIDS))
run_tissue <- function(tis) {
  cc <- grep(sprintf("^%s\\d+_WTS$", tis), names(fpkm), value = TRUE); es <- gsub("\\D", "", cc)
  keep <- es %in% PAIRIDS; cc <- cc[keep]; es <- es[keep]
  dec <- fread(file.path(V, sprintf("decon_v2/CPGEA_%s_cell_fraction.tsv", if (tis == "T") "tumor" else "normal")))
  dec <- dec[, .(sid = gsub("\\D", "", get(names(dec)[1])), epi = Epithelial_cells)]
  cv <- merge(data.table(sid = es), dec, by = "sid", all.x = TRUE, sort = FALSE)
  if (tis == "T") cv <- merge(cv, pur, by = "sid", all.x = TRUE, sort = FALSE)
  cv <- cv[match(es, sid)]; ok <- complete.cases(cv); cc <- cc[ok]; es <- es[ok]; cv <- cv[ok]; n <- length(es)
  E <- log2(as.matrix(fpkm[, ..cc]) + 1); rownames(E) <- fpkm$eg
  ex <- E[rowMeans(E > 0.5) > 0.5, , drop = FALSE]
  pcs <- prcomp(t(ex - rowMeans(ex)), rank. = 5)$x
  base <- cbind(1, pcs, epi = cv$epi, if (tis == "T") cbind(purity = cv$purity) else NULL)
  tfv <- map[, .(sym)][, med := sapply(sym, function(s) { r <- fpkm[gene_name == s, ..cc]; if (!nrow(r)) NA_real_ else max(apply(as.matrix(r), 1, median)) })][is.finite(med) & med >= 1]
  TF <- sapply(tfv$sym, function(s) { r <- as.matrix(fpkm[gene_name == s, ..cc]); r <- r[which.max(apply(r, 1, median)), ]; as.numeric(scale(log2(r + 0.01))) })
  b2 <- bim[, .(idx = idx[1], a1 = a1[1], a2 = a2[1], nd = .N), by = .(chr = as.character(chr), pos)][nd == 1]
  q <- GL[eg %in% rownames(E)]
  q <- merge(q, b2, by = c("chr", "pos"))[(a1 == ref & a2 == alt) | (a1 == alt & a2 == ref)][, .(eg, chr, Pos = pos, idx)]
  GTm <- sapply(q$idx, readsnp); rownames(GTm) <- sids; GTm <- GTm[es, , drop = FALSE]
  msg("[%s] samples %d ; expressed TFs %d ; lead eQTLs with genotype %d", tis, n, ncol(TF), nrow(q))
  Y <- apply(E[q$eg, , drop = FALSE], 1, int)
  one <- function(j, TFm, gxe) {
    G <- GTm[, j]; G[is.na(G)] <- mean(G, na.rm = TRUE); if (var(G) == 0) return(NULL)
    C <- cbind(base, G = G, if (gxe) cbind(Gepi = G * cv$epi) else NULL); Q <- qr.Q(qr(C)); pr <- function(M) M - Q %*% crossprod(Q, M)
    y <- Y[, j]; ry <- as.numeric(pr(y)); sg <- 1
    X1 <- pr(TFm); X2 <- pr(G * TFm)
    a11 <- colSums(X1^2); a22 <- colSums(X2^2); a12 <- colSums(X1 * X2); b1 <- colSums(X1 * ry); b2 <- colSums(X2 * ry)
    dt <- a11 * a22 - a12^2; be2 <- (a11 * b2 - a12 * b1) / dt; be1 <- (a22 * b1 - a12 * b2) / dt
    rss <- sum(ry^2) - (be1 * b1 + be2 * b2); df <- n - ncol(C) - 2
    t <- be2 / sqrt(rss / df * a11 / dt); t[!is.finite(t)] <- NA; t * sg }
  summ <- function(TFm, gxe, keep = FALSE) {
    tt <- lapply(seq_len(nrow(q)), function(j) one(j, TFm, gxe)); v <- unlist(tt)
    out <- c(mean_t = mean(v, na.rm = TRUE), lambda = median(v^2, na.rm = TRUE) / 0.455, frac3 = mean(abs(v) > 3, na.rm = TRUE), n = sum(is.finite(v)))
    if (keep) attr(out, "gene") <- sapply(tt, function(x) if (is.null(x)) NA else mean(x, na.rm = TRUE)); out }
  lam <- function(TFm) { v <- unlist(lapply(seq_len(nrow(q)), function(j) one(j, TFm, TRUE))); c(lambda = median(v^2, na.rm = TRUE) / 0.455, frac3 = mean(abs(v) > 3, na.rm = TRUE)) }
  allm <- apply(E, 1, median); cand <- setdiff(names(allm)[allm >= log2(2) & !duplicated(names(allm))], fpkm[gene_name %in% colnames(TF), eg])
  tfmed <- apply(E[fpkm[match(colnames(TF), gene_name), eg], , drop = FALSE], 1, median)
  ob <- lam(TF)
  RG <- do.call(rbind, mclapply(seq_len(NDRAW), function(k) { set.seed(500 + k)
    br <- cut(allm[cand], quantile(tfmed, 0:5 / 5, na.rm = TRUE), include.lowest = TRUE); need <- table(cut(tfmed, quantile(tfmed, 0:5 / 5, na.rm = TRUE), include.lowest = TRUE))
    pick <- unlist(lapply(names(need), function(b) sample(cand[which(br == b)], min(need[[b]], sum(br == b, na.rm = TRUE)))))
    M <- t(scale(t(E[pick, , drop = FALSE]))); lam(t(M)) }, mc.cores = NC))
  msg("[%s] M1b lambda: TFs %.3f (|t|>3 %.4f) ; expression-matched random non-TF genes %.3f +- %.3f (range %.3f-%.3f; |t|>3 %.4f) ; TFs above random in %d/%d draws",
      tis, ob["lambda"], ob["frac3"], mean(RG[, 1]), sd(RG[, 1]), min(RG[, 1]), max(RG[, 1]), mean(RG[, 2]), sum(ob["lambda"] > RG[, 1]), NDRAW)
  fwrite(data.table(tissue = tis, draw = seq_len(NDRAW), lambda = RG[, 1], frac3 = RG[, 2], tf_lambda = ob["lambda"], tf_frac3 = ob["frac3"]), file.path(O, sprintf("random_nonTF_draws_%s.csv", tis)))
  msg("[%s] empirical one-sided P (TF lambda <= random) = %.4f over %d draws ; TF - random mean = %+.4f", tis, (1 + sum(RG[, 1] >= ob["lambda"])) / (NDRAW + 1), NDRAW, ob["lambda"] - mean(RG[, 1]))
  key <- file.path(O, sprintf("lead_keys_%s.tsv", tis)); fwrite(unique(q[, .(k = paste0("chr", chr, ":", Pos))]), key, col.names = FALSE)
  af <- file.path(O, sprintf("adastra_lead_%s.tsv", tis))
  system(sprintf(paste0("awk -F'\\t' -v OFS='\\t' 'NR==FNR{k[$1]=1;next} FNR==1{for(i=1;i<=NF;i++)c[$i]=i; n=split(FILENAME,a,\"/\"); tf=a[n]; sub(\"_HUMAN.tsv\",\"\",tf); next}",
    " (($1\":\"$3) in k){fr=$(c[\"fdrp_bh_ref\"]); fa=$(c[\"fdrp_bh_alt\"]); if((fr!=\"\" && fr<0.05) || (fa!=\"\" && fa<0.05)) print tf,$1,$3,$5,$6,fr,fa}' %s %s/*_HUMAN.tsv > %s"), key, AD, af))
  A <- fread(af, header = FALSE, col.names = c("mn", "chr", "pos", "ref", "alt", "fdr_ref", "fdr_alt"))
  A <- merge(A, map[, .(mn, tf = sym)], by = "mn")[tf %in% colnames(TF)]
  A[, `:=`(chr = sub("^chr", "", chr), pref_ref = fifelse(is.na(fdr_alt), TRUE, fifelse(is.na(fdr_ref), FALSE, fdr_ref < fdr_alt)))]
  qq <- q[, .(j = .I, chr, Pos, idx)]; qq <- merge(qq, bim[, .(idx, a1, a2)], by = "idx")
  A <- merge(A, qq, by.x = c("chr", "pos"), by.y = c("chr", "Pos"))
  A <- A[(ref == a1 & alt == a2) | (ref == a2 & alt == a1)]
  A[, pref_is_a1 := fifelse(pref_ref, ref == a1, alt == a1)]
  if (nrow(A) >= 5) {
    tt <- lapply(sort(unique(A$j)), function(j) { t <- one(j, TF, TRUE); data.table(j = j, tf = colnames(TF), t = t) }); tt <- rbindlist(tt)
    A <- merge(A[, .(j, tf, pref_is_a1)], tt, by = c("j", "tf"))[is.finite(t)]
    A[, ts := t * fifelse(pref_is_a1, 1, -1)]
    other <- tt[!paste(j, tf) %in% A[, paste(j, tf)] & is.finite(t)]
    sh <- replicate(2000, { s <- sample(A$ts * sample(c(-1, 1), nrow(A), replace = TRUE)); mean(s) })
    msg("[%s] ASB at lead eQTL: %d gene x TF pairs (%d genes, %d TFs) ; ASB-signed interaction t mean %+.3f (sign-flip null P=%.3f) ; mean t^2 %.3f vs other TFs at the same genes %.3f (Wilcoxon P=%.3g)",
        tis, nrow(A), uniqueN(A$j), uniqueN(A$tf), mean(A$ts), (1 + sum(abs(sh) >= abs(mean(A$ts)))) / 2001, mean(A$t^2), mean(other$t^2), wilcox.test(A$t^2, other$t^2)$p.value)
    fwrite(A, file.path(O, sprintf("M1_asb_lead_pairs_%s.csv", tis)))
  } else msg("[%s] ASB at lead eQTL: only %d pairs", tis, nrow(A))
  invisible(NULL) }
invisible(lapply(c("T", "N"), run_tissue)); msg("DONE")
