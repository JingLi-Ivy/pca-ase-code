#!/usr/bin/env Rscript
# ============================================================================
# 09_nested_risk_model.R
# Nested prognostic risk models: clinical (M1) vs + expression-dysregulation burden (M2) vs + ASE burden (M3)
#
# Question: is ASE merely a proxy for expression? Tests whether the global ASE burden adds
# prognostic value beyond clinical factors plus the expression-dysregulation burden of the same gene set.
#
# Endpoint: BCR-free survival (BCR=YES is the event; Dead/NO are censored; time = Month to surgery/BCR)
# Evaluation: likelihood-ratio test / C-index (bootstrap optimism-corrected) / delta-C bootstrap CI /
#       time-dependent AUC (timeROC) / NRI&IDI (survIDINRI) / calibration curve / DCA
#
# Raw data are read-only; results are written to nested_risk_model_results/
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survivalROC)
  library(survIDINRI)
  library(dcurves)
})

set.seed(20260612)

CLIN_FILE <- "/path/to/data_root/CPGEA/CPGEA_web/CPGEA_Clinical_Info.tsv"
FPKM_FILE <- "/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls"
ASE_FILE  <- "/path/to/data_root/ASE/all_T_samples_SNPs_with_ASE_annotation.csv"
PCAASE_FILE <- "/path/to/data_root/ASE/code/PCa_ASE_genes.csv"
OUT_DIR   <- "/path/to/data_root/ASE/nested_risk_model_results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

strip_ver <- function(x) sub("\\..*$", "", x)
log_msg <- function(...) cat(sprintf(...), "\n")

log_msg("=== 1. Reading clinical data ===")
clin_raw <- fread(CLIN_FILE)
clin <- clin_raw[, .(
  sample_id = as.character(`Sequencing ID`),
  age       = as.numeric(Age),
  psa       = as.numeric(PSA),
  gleason   = as.numeric(GS),
  path_t_raw= as.character(`Pathology T Stage`),
  bcr       = as.character(BCR),
  months_str= as.character(`Month to surgery/BCR`)
)]
clin[, path_t := as.numeric(sub("^[^0-9]*([0-9]).*$", "\\1", path_t_raw))]
clin[, bcr_event := fifelse(toupper(bcr) == "YES", 1L,
                     fifelse(bcr %in% c("NA", "", NA), NA_integer_, 0L))]
clin[, surv_time := as.numeric(gsub("[^0-9.]", "", months_str))]
log_msg("Clinical records: %d rows", nrow(clin))

log_msg("=== 2. PCa-ASE gene set ===")
pca_ase <- fread(PCAASE_FILE)
ase_geneset <- unique(strip_ver(pca_ase$gene_id))
log_msg("Number of PCa-ASE genes: %d", length(ase_geneset))

log_msg("=== 3. Computing ASE burden (reading the all-sample SNP annotation, 3 columns only) ===")
ase <- fread(ASE_FILE, select = c("sample_id", "gene_id", "ase_status"))
ase[, gene_id := strip_ver(gene_id)]
ase <- ase[gene_id %in% ase_geneset]
ase[, sid := sub("^[TN]", "", sample_id)]
ase[, is_sig := as.integer(ase_status == "Significant_ASE")]
gene_lv <- ase[, .(sig = as.integer(any(is_sig == 1L))), by = .(sid, gene_id)]
ase_burden <- gene_lv[, .(
  n_tested_ase = .N,
  n_sig_ase    = sum(sig),
  ase_burden   = sum(sig) / .N
), by = sid]
setnames(ase_burden, "sid", "sample_id")
log_msg("Tumor samples with ASE data: %d", nrow(ase_burden))
rm(ase, gene_lv); gc()

log_msg("=== 4. Computing expression-dysregulation burden (paired T/N FPKM) ===")
fpkm <- fread(FPKM_FILE)
fpkm[, gid := strip_ver(gene_id)]
fpkm <- fpkm[gid %in% ase_geneset]
id_cols <- intersect(c("gene_id", "gene_name", "gid"), names(fpkm))
val_cols <- setdiff(names(fpkm), id_cols)
fl <- melt(fpkm, id.vars = "gid", measure.vars = val_cols,
           variable.name = "col", value.name = "fpkm")
fl[, tissue := substr(as.character(col), 1, 1)]
fl[, sid := sub("^[TN]", "", sub("_WTS$", "", as.character(col)))]
fl <- fl[tissue %in% c("T", "N")]
fw <- dcast(fl, sid + gid ~ tissue, value.var = "fpkm", fun.aggregate = mean)
fw <- fw[!is.na(T) & !is.na(N)]
fw[, l2fc := log2((T + 1) / (N + 1))]
fw[, dysreg := as.integer(abs(l2fc) > 1)]
expr_burden <- fw[, .(
  n_tested_expr = .N,
  n_dysreg      = sum(dysreg),
  expr_burden   = sum(dysreg) / .N
), by = sid]
setnames(expr_burden, "sid", "sample_id")
log_msg("Tumor samples with paired expression: %d", nrow(expr_burden))
rm(fpkm, fl, fw); gc()

log_msg("=== 5. Merging ===")
dt <- Reduce(function(a, b) merge(a, b, by = "sample_id"),
             list(clin, ase_burden, expr_burden))
dt <- dt[!is.na(surv_time) & surv_time > 0 & !is.na(bcr_event) &
         !is.na(age) & !is.na(psa) & !is.na(gleason) & !is.na(path_t)]
dt[, ase_burden_z  := scale(ase_burden)[, 1]]
dt[, expr_burden_z := scale(expr_burden)[, 1]]
log_msg("Final modelling samples: n=%d, BCR events=%d (EPV ~ %.1f)",
        nrow(dt), sum(dt$bcr_event), sum(dt$bcr_event) / 6)
fwrite(dt, file.path(OUT_DIR, "nested_model_data.csv"))

if (sum(dt$bcr_event) < 15)
  log_msg("!! Warning: fewer than 15 events; the multivariable model has limited power, interpret with caution")

log_msg("=== 6. Nested Cox models ===")
S <- Surv(dt$surv_time, dt$bcr_event)
f1 <- as.formula("S ~ age + psa + gleason + path_t")
f2 <- update(f1, . ~ . + expr_burden_z)
f3 <- update(f2, . ~ . + ase_burden_z)
m1 <- coxph(f1, data = dt); m2 <- coxph(f2, data = dt); m3 <- coxph(f3, data = dt)

cat("\n--- M1 clinical ---\n");                print(summary(m1)$coefficients)
cat("\n--- M2 + expression-dysregulation burden ---\n");     print(summary(m2)$coefficients)
cat("\n--- M3 +ASE burden ---\n");          print(summary(m3)$coefficients)

lrt <- function(a, b) {
  stat <- 2 * (as.numeric(logLik(b)) - as.numeric(logLik(a)))
  ddf  <- length(coef(b)) - length(coef(a))
  data.table(stat = stat, df = ddf, p = pchisq(stat, ddf, lower.tail = FALSE))
}
lr_tab <- rbindlist(list(
  cbind(comparison = "M2 vs M1 (+expr)", lrt(m1, m2)),
  cbind(comparison = "M3 vs M2 (+ASE)",  lrt(m2, m3)),
  cbind(comparison = "M3 vs M1 (+both)", lrt(m1, m3))
))
cat("\n--- Likelihood-ratio test ---\n"); print(lr_tab)

log_msg("=== 7. C-index (apparent + optimism-corrected) ===")
cstat <- function(fit, newdata = NULL) {
  if (is.null(newdata)) as.numeric(concordance(fit)$concordance)
  else as.numeric(concordance(fit, newdata = newdata)$concordance)
}
c_app <- c(cstat(m1), cstat(m2), cstat(m3))

log_msg("bootstrap (optimism + ΔC) ...")
B <- 1000
opt  <- matrix(NA, B, 3)
dC   <- matrix(NA, B, 3, dimnames = list(NULL, c("M2-M1", "M3-M2", "M3-M1")))
for (b in seq_len(B)) {
  idx <- sample.int(nrow(dt), replace = TRUE)
  db <- dt[idx]
  fb <- lapply(list(f1, f2, f3), function(ff)
          tryCatch(coxph(ff, data = db), error = function(e) NULL))
  if (any(vapply(fb, is.null, logical(1)))) next
  c_boot <- vapply(fb, cstat, numeric(1))
  c_test <- vapply(fb, function(f) tryCatch(cstat(f, newdata = dt),
                    error = function(e) NA_real_), numeric(1))
  opt[b, ] <- c_boot - c_test
  dC[b, ]  <- c(c_boot[2] - c_boot[1], c_boot[3] - c_boot[2], c_boot[3] - c_boot[1])
}
c_corrected <- c_app - colMeans(opt, na.rm = TRUE)
c_tab <- data.table(
  model       = c("M1_clinical", "M2_+expr", "M3_+ASE"),
  c_apparent  = c_app,
  optimism    = colMeans(opt, na.rm = TRUE),
  c_corrected = c_corrected
)
cat("\n--- C-index (B=1000 optimism-corrected) ---\n"); print(c_tab)

dC_tab <- data.table(
  comparison = colnames(dC),
  dC_mean = colMeans(dC, na.rm = TRUE),
  ci_lo   = apply(dC, 2, quantile, 0.025, na.rm = TRUE),
  ci_hi   = apply(dC, 2, quantile, 0.975, na.rm = TRUE),
  p_oneside = apply(dC, 2, function(x) mean(x <= 0, na.rm = TRUE))
)
cat("\n--- ΔC bootstrap (B=1000) ---\n"); print(dC_tab)

log_msg("=== 8. Time-dependent AUC ===")
lp <- data.table(M1 = predict(m1, type = "lp"),
                 M2 = predict(m2, type = "lp"),
                 M3 = predict(m3, type = "lp"))
times <- c(12, 36, 60)
times <- times[times < max(dt$surv_time)]
auc_list <- lapply(names(lp), function(nm) {
  aucs <- vapply(times, function(tt)
    survivalROC(Stime = dt$surv_time, status = dt$bcr_event,
                marker = lp[[nm]], predict.time = tt, method = "KM")$AUC,
    numeric(1))
  data.table(model = nm, time = times, AUC = aucs)
})
auc_tab <- rbindlist(auc_list)
cat("\n--- Time-dependent AUC ---\n"); print(dcast(auc_tab, time ~ model, value.var = "AUC"))

log_msg("=== 9. NRI & IDI ===")
t0 <- if (60 < quantile(dt$surv_time, 0.9)) 60 else round(median(dt$surv_time))
outcome <- as.matrix(dt[, .(surv_time, bcr_event)])
X1 <- model.matrix(f1, dt)[, -1, drop = FALSE]
X2 <- model.matrix(f2, dt)[, -1, drop = FALSE]
X3 <- model.matrix(f3, dt)[, -1, drop = FALSE]
nri_idi <- function(xold, xnew, label) {
  r <- IDI.INF(outcome, xold, xnew, t0, npert = 200)
  data.table(comparison = label, t0 = t0,
             IDI = r$m1[1], IDI_lo = r$m1[2], IDI_hi = r$m1[3], IDI_p = r$m1[4],
             NRI = r$m2[1], NRI_lo = r$m2[2], NRI_hi = r$m2[3], NRI_p = r$m2[4])
}
ni_tab <- tryCatch(rbindlist(list(
  nri_idi(X1, X2, "M2 vs M1 (+expr)"),
  nri_idi(X2, X3, "M3 vs M2 (+ASE)"),
  nri_idi(X1, X3, "M3 vs M1 (+both)")
)), error = function(e) { log_msg("NRI/IDI failed: %s", conditionMessage(e)); NULL })
if (!is.null(ni_tab)) { cat("\n--- NRI & IDI @t0=", t0, "mo ---\n"); print(ni_tab) }

log_msg("=== 10. Calibration curve + DCA ===")
risk_at <- function(fit, u) {
  sf <- survfit(fit, newdata = dt)
  s_u <- summary(sf, times = u)$surv
  1 - as.numeric(s_u)
}
u <- t0
dt[, `:=`(risk_M1 = risk_at(m1, u), risk_M2 = risk_at(m2, u), risk_M3 = risk_at(m3, u))]

pdf(file.path(OUT_DIR, "Fig_nested_risk_model.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2))
bp <- barplot(rbind(c_tab$c_apparent, c_tab$c_corrected), beside = TRUE,
        names.arg = c_tab$model, ylim = c(0.4, 1), col = c("grey70", "steelblue"),
        legend.text = c("Apparent", "Optimism-corrected"),
        main = "(a) Discrimination (C-index)", ylab = "C-index")
abline(h = 0.5, lty = 2)
plot(NA, xlim = range(times), ylim = c(0.4, 1), xlab = "Month", ylab = "AUC(t)",
     main = "(b) Time-dependent AUC")
cols <- c(M1 = "grey50", M2 = "darkorange", M3 = "firebrick")
for (nm in names(lp)) {
  sub <- auc_tab[model == nm]; lines(sub$time, sub$AUC, col = cols[nm], lwd = 2, type = "b", pch = 19)
}
abline(h = 0.5, lty = 2); legend("bottomright", names(cols), col = cols, lwd = 2, pch = 19)
K <- 4
grp <- cut(dt$risk_M3, breaks = quantile(dt$risk_M3, probs = seq(0, 1, length.out = K + 1)),
           include.lowest = TRUE, labels = FALSE)
cal_tab <- rbindlist(lapply(sort(unique(grp)), function(g) {
  sub <- dt[grp == g]
  km <- survfit(Surv(surv_time, bcr_event) ~ 1, data = sub)
  ss <- summary(km, times = u, extend = TRUE)
  data.table(group = g, pred = mean(sub$risk_M3),
             obs = 1 - ss$surv, obs_lo = 1 - ss$upper, obs_hi = 1 - ss$lower)
}))
plot(cal_tab$pred, cal_tab$obs, xlim = c(0, 1), ylim = c(0, 1), pch = 19, col = "firebrick",
     xlab = "Predicted risk", ylab = "Observed (KM) risk",
     main = sprintf("(c) Calibration M3 @ %d mo", u))
arrows(cal_tab$pred, cal_tab$obs_lo, cal_tab$pred, cal_tab$obs_hi,
       angle = 90, code = 3, length = 0.05, col = "firebrick")
abline(0, 1, lty = 2)
dev.off()
fwrite(cal_tab, file.path(OUT_DIR, "calibration_M3.csv"))

dca_dt <- dt[, .(surv_time, bcr_event, risk_M1, risk_M2, risk_M3)]
dca_res <- tryCatch(
  dca(Surv(surv_time, bcr_event) ~ risk_M1 + risk_M2 + risk_M3,
      data = dca_dt, time = u),
  error = function(e) { log_msg("DCA failed: %s", conditionMessage(e)); NULL })
if (!is.null(dca_res)) {
  ggplot2::ggsave(file.path(OUT_DIR, "Fig_nested_DCA.pdf"),
                  plot(dca_res, smooth = TRUE) + ggplot2::ggtitle(sprintf("DCA @ %d mo", u)),
                  width = 7, height = 5)
}

fwrite(lr_tab, file.path(OUT_DIR, "LR_tests.csv"))
fwrite(c_tab,  file.path(OUT_DIR, "Cindex.csv"))
fwrite(dC_tab, file.path(OUT_DIR, "deltaC_bootstrap.csv"))
fwrite(auc_tab,file.path(OUT_DIR, "timeAUC.csv"))
if (!is.null(ni_tab)) fwrite(ni_tab, file.path(OUT_DIR, "NRI_IDI.csv"))
coef_out <- rbindlist(lapply(list(M1=m1,M2=m2,M3=m3), function(f){
  s <- summary(f)$coefficients; ci <- summary(f)$conf.int
  data.table(term=rownames(s), HR=ci[,1], lo=ci[,3], hi=ci[,4], p=s[,5])
}), idcol = "model")
fwrite(coef_out, file.path(OUT_DIR, "model_coefficients.csv"))

log_msg("=== Done. Results written to %s ===", OUT_DIR)
sink(file.path(OUT_DIR, "session_info.txt")); print(sessionInfo()); sink()
