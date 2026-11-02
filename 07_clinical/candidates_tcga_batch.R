#!/usr/bin/env Rscript
# TCGA-PRAD validation for all 6 candidates
suppressPackageStartupMessages({library(data.table); library(survival)})

load('/path/to/data_root/ASE/code/fig4_results/tcga_validation/TCGA_validation_workspace.RData')

clin_dt <- as.data.table(clinical_matched, keep.rownames = "patient_id")

clin_bcr <- clin_dt[, .(patient_id)]
clin_bcr[, bcr_days_val := as.numeric(clin_dt$patient.days_to_first_biochemical_recurrence)]
clin_bcr[, fu_days := as.numeric(clin_dt$patient.days_to_last_followup)]
clin_bcr[, bcr_event := as.integer(!is.na(bcr_days_val) & bcr_days_val > 0)]
clin_bcr[, time_days := ifelse(bcr_event == 1, bcr_days_val, fu_days)]
clin_bcr[, surv_months := time_days / 30.44]
clin_bcr <- clin_bcr[!is.na(surv_months) & surv_months > 0 & !is.na(bcr_event)]
cat(sprintf("TCGA-PRAD: %d patients, %d BCR events\n\n", nrow(clin_bcr), sum(clin_bcr$bcr_event)))

suppressPackageStartupMessages({library(curatedTCGAData); library(MultiAssayExperiment)})
prad <- curatedTCGAData("PRAD", "RNASeq2Gene", version = "2.0.1", dry.run = FALSE)
expr <- assay(prad[["PRAD_RNASeq2Gene-20160128"]])

gs_cols <- grep("gleason", names(clin_dt), value = TRUE, ignore.case = TRUE)
cat(sprintf("Gleason columns found: %s\n", paste(gs_cols, collapse=", ")))
if (length(gs_cols) > 0) {
  gs_col <- gs_cols[1]
  gs_vals <- clin_dt[[gs_col]][match(clin_bcr$patient_id, clin_dt$patient_id)]
  clin_bcr[, gleason_score := as.numeric(gs_vals)]
  clin_bcr[, gleason_high := as.integer(gleason_score >= 8)]
  cat(sprintf("Gleason: mean=%.1f, high(>=8): %d/%d\n",
              mean(clin_bcr$gleason_score, na.rm=T), sum(clin_bcr$gleason_high==1, na.rm=T), nrow(clin_bcr)))
} else {
  clin_bcr[, gleason_high := NA_integer_]
}

candidates <- c("CELF2", "PPP3CA", "HSPH1", "ZBTB16", "ITPR1", "SERINC5")

avail <- candidates[candidates %in% rownames(expr)]
missing <- candidates[!candidates %in% rownames(expr)]
if (length(missing) > 0) cat(sprintf("WARNING: Not in TCGA expression: %s\n", paste(missing, collapse=", ")))
cat(sprintf("Available: %s\n\n", paste(avail, collapse=", ")))

for (gene in avail) {
  cat(sprintf("\n%s\n", paste(rep("-", 50), collapse="")))
  cat(sprintf("--- TCGA Validation: %s ---\n", gene))
  cat(sprintf("%s\n", paste(rep("-", 50), collapse="")))

  gene_expr <- data.table(sample_id = colnames(expr), value = as.numeric(expr[gene, ]))
  gene_expr[, patient_id := sapply(strsplit(sample_id, "-"), function(x) paste(x[1:3], collapse="-"))]
  gene_expr <- gene_expr[grepl("-01[A-Z]", sample_id)]
  gene_expr <- gene_expr[, .SD[1], by = patient_id]

  surv_df <- merge(gene_expr, clin_bcr, by = "patient_id")
  surv_df[, expr_log := log2(value + 1)]

  cat(sprintf("  Matched: %d patients, %d events\n", nrow(surv_df), sum(surv_df$bcr_event)))
  cat(sprintf("  Expression: mean=%.2f, median=%.2f (log2)\n", mean(surv_df$expr_log), median(surv_df$expr_log)))

  cox1 <- coxph(Surv(surv_months, bcr_event) ~ expr_log, data = surv_df)
  s1 <- summary(cox1)
  cat(sprintf("  Cox univariate: HR=%.3f (%.3f-%.3f), p=%.4f\n",
              s1$conf.int[1,1], s1$conf.int[1,3], s1$conf.int[1,4], s1$coefficients[1,5]))

  surv_multi <- surv_df[!is.na(gleason_high)]
  if (nrow(surv_multi) > 30) {
    cox2 <- coxph(Surv(surv_months, bcr_event) ~ expr_log + gleason_high, data = surv_multi)
    s2 <- summary(cox2)
    cat(sprintf("  Cox multivariate (adj Gleason): HR=%.3f, p=%.4f\n",
                s2$conf.int[1,1], s2$coefficients[1,5]))
  }

  med <- median(surv_df$expr_log)
  surv_df[, group := ifelse(expr_log >= med, "High", "Low")]
  sdiff <- survdiff(Surv(surv_months, bcr_event) ~ group, data = surv_df)
  lr_p <- 1 - pchisq(sdiff$chisq, 1)

  ev_low <- sum(surv_df$group == "Low" & surv_df$bcr_event == 1)
  ev_high <- sum(surv_df$group == "High" & surv_df$bcr_event == 1)
  n_low <- sum(surv_df$group == "Low")
  n_high <- sum(surv_df$group == "High")
  cat(sprintf("  Log-rank: p=%.4f\n", lr_p))
  cat(sprintf("    Low: n=%d, events=%d (%.1f%%)\n", n_low, ev_low, 100*ev_low/n_low))
  cat(sprintf("    High: n=%d, events=%d (%.1f%%)\n", n_high, ev_high, 100*ev_high/n_high))

  cat(sprintf("  Direction (HR<1 = low expr worse): %s\n",
              ifelse(s1$conf.int[1,1] < 1, "LOW worse (consistent)", "HIGH worse")))
  cat("\n")
}

cat("\n=== TCGA BATCH VALIDATION COMPLETE ===\n")
