#!/usr/bin/env Rscript
# A11: external validation of the locked KIF18B expression -> BCR model (analysis plan pre-specified on 2026-09-26)
suppressPackageStartupMessages({library(data.table); library(survival)})
setwd("/path/to/data_root/ASE/code/v8_round2/external")
gs <- function(x) { x <- as.character(x); s <- suppressWarnings(as.numeric(x)); i <- grepl("\\+", x)
  s[i] <- sapply(strsplit(x[i], "+", fixed = TRUE), function(v) sum(as.numeric(v))); s }
ser <- function(f) {
  L <- readLines(gzfile(f)); ch <- L[startsWith(L, "!Sample_characteristics_ch1")]
  ids <- strsplit(L[startsWith(L, "!Sample_geo_accession")], "\t")[[1]][-1]
  M <- lapply(ch, function(l) gsub('"', "", strsplit(l, "\t")[[1]][-1]))
  key <- sapply(M, function(v) sub(":.*", "", v[1]))
  d <- as.data.table(setNames(lapply(M, function(v) sub("^[^:]*: ?", "", v)), make.names(key, unique = TRUE)))
  d[, gsm := gsub('"', "", ids)]; list(d = d, L = L) }
C <- list()

e <- fread("prad_mskcc_KIF18B.csv"); c <- fread("prad_mskcc_clinical.csv")
d <- merge(e, c, by = c("sampleId","patientId"))[SAMPLE_TYPE == "Primary" & DFS_STATUS != ""]
C$MSKCC_Taylor2010 <- d[, .(id = patientId, expr = KIF18B, time = as.numeric(DFS_MONTHS), event = as.integer(DFS_STATUS == "1:Recurred"),
                             gleason = as.numeric(GLEASON_SCORE), psa = NA_real_)]
e <- fread("prostate_dkfz_2018_KIF18B.csv"); c <- fread("prostate_dkfz_2018_clinical.csv")
d <- merge(e, c, by = c("sampleId","patientId"))[INITIAL_TREATMENT == "RadP" & !is.na(BCR_STATUS)]
C$DKFZ_2018 <- d[, .(expr = mean(log2(KIF18B + 0.01)), time = as.numeric(TIME_FROM_SURGERY_TO_BCR_LASTFU[1]), event = as.integer(BCR_STATUS[1]),
                     gleason = gs(GLEASON_SCORE[1]), psa = as.numeric(PREOPERATIVE_PSA[1])), by = .(id = patientId)]
s <- ser("GSE54460_series_matrix.txt.gz")$d
f <- fread("GSE54460_FPKM-genes-TopHat2-106samples-12-4-13.txt.gz", header = FALSE)
ext <- unlist(f[2, -(1:2)]); kv <- as.numeric(unlist(f[V2 == "KIF18B", -(1:2)]))
k <- data.table(ext = ext, expr = log2(kv + 0.01))
s[, ext := external.id]; s <- merge(s, k, by = "ext")
C$GSE54460_Long2014 <- s[, .(id = ext, expr, event = as.integer(bcr),
    time = fifelse(bcr == "1", suppressWarnings(as.numeric(months.to.bcr)), suppressWarnings(as.numeric(months.total.f.u))),
    gleason = as.numeric(tgleason), psa = suppressWarnings(as.numeric(prepsa)))]
if (file.exists("GSE116918_KIF18B_probes.txt")) {
  pr <- readLines("GSE116918_KIF18B_probes.txt"); S <- ser("GSE116918_series_matrix.txt.gz")
  L <- S$L; b <- which(L == "!series_matrix_table_begin"); t <- fread(text = L[(b + 1):(which(L == "!series_matrix_table_end") - 1)])
  x <- t[get(names(t)[1]) %in% pr]; ev <- colMeans(as.matrix(x[, -1]))
  s <- S$d; s[, expr := ev[gsm]]
  C$GSE116918_Jain2018 <- s[, .(id = gsm, expr, event = as.integer(bcr.event..1.yes..0.no.), time = as.numeric(follow.up.time..bcr..months.),
                                gleason = as.numeric(gleason.grade), psa = suppressWarnings(as.numeric(psa..ng.ml.)))]
  cat("GSE116918 KIF18B probes:", length(pr), "\n") }

cidx <- function(fit) summary(fit)$concordance[1]
out <- rbindlist(lapply(names(C), function(nm) {
  d <- C[[nm]][is.finite(expr) & is.finite(time) & time > 0 & !is.na(event)]
  d[, z := as.numeric(scale(expr))][, hi := as.integer(expr > median(expr))][, g8 := as.integer(gleason >= 8)]
  u <- coxph(Surv(time, event) ~ z, d); su <- summary(u)$coef
  lr <- survdiff(Surv(time, event) ~ hi, d); plr <- 1 - pchisq(lr$chisq, 1)
  hm <- summary(coxph(Surv(time, event) ~ hi, d))$conf.int
  usepsa <- sum(is.finite(d$psa)) >= 0.8 * nrow(d)
  dd <- if (usepsa) d[is.finite(psa) & psa > 0 & is.finite(gleason)] else d[is.finite(gleason)]
  f0 <- if (usepsa) Surv(time, event) ~ g8 + log(psa) else Surv(time, event) ~ g8
  f1 <- update(f0, . ~ . + z)
  a <- summary(coxph(f1, dd))$coef["z", ]
  c0 <- cidx(coxph(f0, dd)); c1 <- cidx(coxph(f1, dd))
  set.seed(1); bd <- replicate(1000, { i <- sample(nrow(dd), replace = TRUE); x <- dd[i]
    tryCatch(cidx(coxph(f1, x)) - cidx(coxph(f0, x)), error = function(e) NA) })
  data.table(cohort = nm, n = nrow(d), events = sum(d$event), zero_expr_frac = round(mean(d$expr <= log2(0.01) + 1e-9), 2),
    HR_perSD = exp(su[1]), lo = exp(su[1] - 1.96 * su[3]), hi = exp(su[1] + 1.96 * su[3]), P = su[5], logHR = su[1], se = su[3],
    adj_covars = if (usepsa) "Gleason>=8 + log PSA" else "Gleason>=8", n_adj = nrow(dd), HR_adj = exp(a[1]), P_adj = a[5],
    C_clin = c0, C_clin_KIF18B = c1, dC_lo = quantile(bd, .025, na.rm = TRUE), dC_hi = quantile(bd, .975, na.rm = TRUE),
    HR_median_split = hm[1, 1], logrank_P = plr) }))
w <- 1 / out$se^2; fe <- sum(w * out$logHR) / sum(w); Q <- sum(w * (out$logHR - fe)^2); k <- nrow(out)
tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))); wr <- 1 / (out$se^2 + tau2)
re <- sum(wr * out$logHR) / sum(wr); sre <- sqrt(1 / sum(wr))
cat("\n===== Cohorts (primary: univariable HR per SD) =====\n")
print(out[, .(cohort, n, events, zero_expr_frac, HR = round(HR_perSD, 2), CI = sprintf("%.2f-%.2f", lo, hi), P = signif(P, 2),
              adj = adj_covars, n_adj, HR_adj = round(HR_adj, 2), P_adj = signif(P_adj, 2),
              C0 = round(C_clin, 3), C1 = round(C_clin_KIF18B, 3), dC_CI = sprintf("%+.3f to %+.3f", dC_lo, dC_hi),
              HR_med = round(HR_median_split, 2), logrank_P = signif(logrank_P, 2))], width = 250)
cat(sprintf("\nRandom-effects pooled: HR per SD = %.2f (95%% CI %.2f-%.2f), P = %.3g ; I2 = %.0f%%, tau2 = %.3f, Q P = %.3g\n",
    exp(re), exp(re - 1.96 * sre), exp(re + 1.96 * sre), 2 * pnorm(-abs(re / sre)), 100 * max(0, (Q - (k - 1)) / Q), tau2, 1 - pchisq(Q, k - 1)))
fwrite(out, "external_validation_KIF18B.csv")
