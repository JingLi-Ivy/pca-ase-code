# 63b (2026-10-13): layout/label pass of 63 for v10 (scientific-visualization checklist); data and statistics unchanged; output pub_figs_v10.
#  S11E '×' in model labels; S14 legend not clipped; S16A points instead of bars on a truncated axis; S16C axis limited to the bulk (one patient at -0.10 stated in the legend).
#!/usr/bin/env Rscript
# 63: publication versions of Supplemental Figs S11 (KIF18B locus and FOXM1 analyses), S12 (siRNA functional assays),
#     S13 (KIF18B prognosis), S14 (single-cell localization), S16 (clinical associations: nested Cox, Delta-ASE).
# Inputs: saved results (paths below); wet-lab values read from wetlab/ASE_raw_data.xlsx. Transwell and colony images are
# not redrawn (representative images placed in Illustrator). Legend numbers -> supp_C_legend_stats.txt.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(jsonlite); library(readxl); library(survival)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
P <- "/path/to/data_root/ASE"; O <- file.path(V, "pub_figs_v10"); cT <- PAL[["tumor"]]; cN <- PAL[["normal"]]; TS <- BASE / .pt
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
fin <- function(g, f, w, h) save_panel(g + plot_annotation(tag_levels = "A") & tags, file.path(O, f), w, h)
XL <- file.path(V, "wetlab/ASE_raw_data.xlsx"); SRC <- file.path(V, "source_data_v9")

ld <- fread(file.path(P, "figures/rs9914422_LD/FigS5_LD_matrix.csv")); L <- melt(ld, id.vars = "snp", variable.name = "snp2", value.name = "r2")
lv <- ld$snp; L[, `:=`(snp = factor(snp, levels = lv), snp2 = factor(snp2, levels = rev(lv)))]
note("S11A: pairwise r2 among %s (CPGEA, n=134); rs9914422-rs8070376 %.2f; with GTEx eQTL rs144461621 %.3f-%.3f", paste(lv, collapse = ", "),
     ld[snp == "rs9914422", rs8070376], min(unlist(ld[snp != "rs144461621", rs144461621])), max(unlist(ld[snp != "rs144461621", rs144461621])))
A11 <- ggplot(L, aes(snp, snp2, fill = r2)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", r2)), size = TS, family = FONT) +
  scale_fill_gradient(low = "white", high = cT, limits = c(0, 1), name = "r\u00b2") + labs(x = NULL, y = NULL) + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
sc <- as.data.table(fromJSON(file.path(V, "locus_fig/motif_allele_scores.json"))); if (ncol(sc) == 5) setnames(sc, c("tf", "id", "G", "C", "d")) else sc <- sc[, 1:5]
sc <- sc[order(-abs(as.numeric(d)))][1:min(.N, 12)][, `:=`(G = as.numeric(G), C = as.numeric(C))][, tf := factor(tf, levels = rev(tf))]
note("S11B: JASPAR motif scores at rs9914422 (G vs C), top %d by |difference|: %s", nrow(sc), sc[, paste(sprintf("%s %.1f/%.1f", tf, G, C), collapse = "; ")])
B11 <- ggplot(sc) + geom_segment(aes(x = G, xend = C, y = tf, yend = tf), colour = "grey60", linewidth = 0.4) + geom_point(aes(G, tf), colour = "grey40", size = 1.2) +
  geom_point(aes(C, tf), colour = cT, size = 1.2) + labs(x = "Motif score, G (grey) and C (red)", y = NULL) + theme_pub()
fp <- fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls"); tc <- grep("^T\\d+_WTS$", names(fp), value = TRUE)
fam <- fread("/path/to/data_root/CPGEA/CPGEA_genotype/combined.134.qc.fam", header = FALSE)
tc <- tc[gsub("\\D", "", tc) %in% gsub("\\D", "", fam$V1)]
gv <- function(g) log2(as.numeric(unlist(fp[gene_name == g, ..tc][1])) + 1)
FX <- rbindlist(lapply(c("KIF18B", "NEK2", "TOP2A"), function(g) data.table(gene = g, foxm1 = gv("FOXM1"), y = gv(g))))
for (g in unique(FX$gene)) { x <- FX[gene == g]; ct <- cor.test(x$foxm1, x$y, method = "spearman", exact = FALSE); note("S11C %s vs FOXM1: n=%d tumors, Spearman rho %.3f, P=%.1e", g, nrow(x), ct$estimate, ct$p.value) }
C11 <- ggplot(FX, aes(foxm1, y)) + geom_point(size = 0.5, colour = PAL[["dark"]], alpha = 0.7) + geom_smooth(method = "lm", formula = y ~ x, colour = cT, linewidth = 0.4, se = FALSE) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) + labs(x = "FOXM1, log2 FPKM", y = "log2 FPKM") + theme_pub() + theme(strip.text = element_text(face = "italic"))
au <- fread(file.path(V, "kif18b_audit/kif18b_patient_audit_all.csv"))[is.finite(FOXM1) & is.finite(KIF18B) & is.finite(rs9914422_dosage)]
au[, geno := factor(rs9914422_dosage, levels = 2:0, labels = c("GG", "GC", "CC"))]
note("S11D: %d tumors with genotype and expression; genotype classes %s; T67 FOXM1 %.1f FPKM", nrow(au), au[, .N, by = geno][, paste(geno, N, collapse = ", ")], au[pid == 67, FOXM1])
D11 <- ggplot(au, aes(FOXM1, KIF18B, colour = geno)) + geom_point(size = 0.6) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4) +
  geom_point(data = au[pid == 67], shape = 21, size = 2.5, colour = "black", stroke = 0.5) +
  scale_colour_manual(values = c("#BDBDBD", "#6BAED6", cT), name = "rs9914422") + labs(x = "FOXM1, FPKM", y = "KIF18B, FPKM") + theme_pub()
sm <- fread(file.path(SRC, "SuppFigS17_sensitivity_models.csv"))
note("S11E: genotype x FOXM1 interaction on KIF18B across models: %s", sm[, paste(sprintf("%s %.3f [%.3f, %.3f] (n=%d)", model, b, lo, hi, n), collapse = "; ")])
sm[, model := gsub(" x ", " \u00d7 ", model)]
E11 <- ggplot(sm, aes(b, factor(model, levels = rev(model)))) + geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.2, linewidth = 0.4) + labs(x = "Interaction coefficient", y = NULL) + theme_pub()
mg <- fread(file.path(V, "kif18b_audit/ase_gleason_matched/matched_genes.csv")); kv <- fread(file.path(V, "kif18b_audit/ase_gleason_matched/kif18b_vs_matched.csv"))[primary == TRUE]
note("S11F: KIF18B ASE-Gleason z %.2f vs %d expression-matched genes (median %.2f); empirical two-sided P=%.3f; percentile %.0f", kv$KIF18B, kv$matched_n, kv$matched_median, kv$P_matched, kv$pct_matched)
F11 <- ggplot(mg, aes(z_gs_dep)) + geom_histogram(bins = 30, fill = "grey70", colour = "white", linewidth = 0.15) + geom_vline(xintercept = kv$KIF18B, colour = cT, linewidth = 0.6) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) + labs(x = "ASE-Gleason z", y = "Matched genes") + theme_pub()
fin((A11 | B11) / C11 / (D11 | E11 | F11) + plot_layout(heights = c(1.2, 0.8, 1)), "SuppFig_S11.pdf", 178, 175)

rx <- function(sh, rng) as.data.frame(read_excel(XL, sheet = sh, range = rng, col_names = FALSE, .name_repair = "minimal"))
grp3 <- function(df, cols, gene) data.table(gene = gene, cond = rep(c("siNC", "si#1", "si#2"), each = 3), y = as.numeric(unlist(df[, cols])))
q1 <- rx("qPCR1", "A1:P6"); KD <- rbind(grp3(q1[4:6, ], 2:4, "KIF18B"), grp3(q1[4:6, ], 8:10, "TOP2A"), grp3(q1[4:6, ], 14:16, "NEK2"))
mi <- rx("Migration", "A1:P5"); MG <- rbind(grp3(mi[3:5, ], 2:4, "KIF18B"), grp3(mi[3:5, ], 8:10, "TOP2A"), grp3(mi[3:5, ], 14:16, "NEK2"))
co <- rx("colony formation", "A1:K5"); CO <- rbind(grp3(co[3:5, ], 2:4, "KIF18B"), grp3(co[3:5, ], 8:10, "TOP2A"))
ck <- as.data.frame(read_excel(XL, sheet = "CCK8", range = "A1:AB17", col_names = FALSE, .name_repair = "minimal"))
ckb <- function(rows, c0, gene) rbindlist(lapply(rows, function(r) data.table(gene = gene, day = as.numeric(ck[r, c0]),
  cond = rep(c("siNC", "si#1", "si#2"), each = 4), y = as.numeric(unlist(ck[r, (c0 + 1):(c0 + 12)])))))
CK <- rbind(ckb(3:6, 1, "KIF18B"), ckb(3:6, 16, "TOP2A"), ckb(14:17, 1, "NEK2"))
for (D in list(list(KD, "knockdown"), list(MG, "migration"), list(CO, "colony"))) note("S12 %s (mean of n=3 wells): %s", D[[2]], D[[1]][, .(m = mean(y)), by = .(gene, cond)][, paste(sprintf("%s %s %.2f", gene, cond, m), collapse = "; ")])
note("S12 CCK-8: n=%s wells per condition and day (legend of v9 S7 said three)", paste(unique(CK[, .N, by = .(gene, day, cond)]$N), collapse = "/"))
lv3 <- c("siNC", "si#1", "si#2"); cols3 <- c(siNC = "grey75", `si#1` = "#4393C3", `si#2` = "#2166AC")
bars <- function(D, ylab) { D[, `:=`(cond = factor(cond, levels = lv3), gene = factor(gene, levels = unique(gene)))]
  M <- D[, .(m = mean(y), s = sd(y)), by = .(gene, cond)]
  ggplot(M, aes(gene, m, fill = cond)) + geom_col(width = 0.7, position = position_dodge(0.75)) +
    geom_errorbar(aes(ymin = m - s, ymax = m + s), width = 0.2, linewidth = 0.3, position = position_dodge(0.75)) +
    scale_fill_manual(values = cols3, name = NULL) + scale_y_continuous(expand = expansion(c(0, 0.05))) + labs(x = NULL, y = ylab) + theme_pub() +
    theme(axis.text.x = element_text(face = "italic")) }
A12 <- bars(KD, "Relative mRNA") + theme(legend.position = "top")
CM <- CK[, .(m = mean(y), s = sd(y)), by = .(gene, day, cond)][, cond := factor(cond, levels = lv3)]
B12 <- ggplot(CM, aes(day, m, colour = cond)) + geom_line(linewidth = 0.4) + geom_errorbar(aes(ymin = m - s, ymax = m + s), width = 0.15, linewidth = 0.3) +
  geom_point(size = 0.8) + facet_wrap(~gene, nrow = 1) + scale_colour_manual(values = cols3, guide = "none") +
  labs(x = "Day", y = "Absorbance, 450 nm") + theme_pub() + theme(strip.text = element_text(face = "italic"))
C12 <- bars(MG, "Relative migrated cells") + theme(legend.position = "none")
D12 <- bars(CO, "Relative colonies") + theme(legend.position = "none")
fin((A12 | B12) / (C12 | D12) + plot_layout(widths = c(1, 1.6)), "SuppFig_S12.pdf", 178, 105)

e <- new.env(); load(file.path(P, "code/fig4_results/tcga_validation/TCGA_validation_workspace.RData"), envir = e)
cm <- e$clinical_matched; fu <- as.numeric(cm[["days_to_last_followup"]]); if (all(is.na(fu))) fu <- as.numeric(cm[["patient.days_to_last_followup"]])
bd <- as.numeric(cm[["patient.days_to_first_biochemical_recurrence"]]); st <- tolower(as.character(cm[["patient.biochemical_recurrence"]]))
ev <- ifelse(st %in% c("yes", "true", "1"), 1, ifelse(st %in% c("no", "false", "0"), 0, NA)); tv <- ifelse(ev == 1, ifelse(is.na(bd), fu, bd), fu)
base <- data.table(id = cm$tcga_full_id, ev = ev, t = tv / 30.44)[!is.na(ev) & !is.na(t) & t > 0 & !is.na(id)]
d <- merge(base, data.table(id = colnames(e$expr_matrix_log2), x = as.numeric(e$expr_matrix_log2["KIF18B", ])), by = "id")
d[, grp := factor(fifelse(x >= median(x), "High", "Low"), levels = c("Low", "High"))]
cx <- summary(coxph(Surv(t, ev) ~ grp, d)); lr <- survdiff(Surv(t, ev) ~ grp, d)
note("S13A TCGA-PRAD: n=%d, %d events; KIF18B high vs low HR %.2f [%.2f, %.2f]; log-rank P=%.1e", nrow(d), sum(d$ev), cx$conf.int[1, 1], cx$conf.int[1, 3], cx$conf.int[1, 4], 1 - pchisq(lr$chisq, 1))
sf <- survfit(Surv(t, ev) ~ grp, d); kmd <- data.table(t = sf$time, s = sf$surv, grp = rep(c("Low", "High"), sf$strata))
kmd <- rbind(data.table(t = 0, s = 1, grp = c("Low", "High")), kmd)
A13 <- ggplot(kmd, aes(t, s, colour = grp)) + geom_step(linewidth = 0.4) + scale_colour_manual(values = c(Low = cN, High = cT), name = "KIF18B") +
  scale_y_continuous(limits = c(0, 1)) + labs(x = "Months", y = "BCR-free survival") + theme_pub() + theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0))
ex <- fread(file.path(V, "external/external_validation_KIF18B.csv"))
ex[, lab := sub("_", " ", cohort)]
note("S13B external cohorts (HR per SD): %s", ex[, paste(sprintf("%s %.2f [%.2f, %.2f] (n=%d, %d events)", lab, HR_perSD, lo, hi, n, events), collapse = "; ")])
B13 <- ggplot(ex, aes(HR_perSD, lab)) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) + geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.25, linewidth = 0.4) +
  scale_x_log10() + labs(x = "HR per SD, log scale", y = NULL) + theme_pub()
bp <- fread(file.path(V, "external/kif18b_beyond_proliferation.csv"))
bp[, `:=`(lab = sub("_.*", "", cohort), sc = fifelse(prolif == "P7", "7-gene score", "CCP score"))]
note("S13C KIF18B HR adjusted for proliferation: %s", bp[, paste(sprintf("%s/%s %.2f [%.2f, %.2f]", lab, sc, HR_K_adjProlif, lo, hi), collapse = "; ")])
C13 <- ggplot(bp, aes(HR_K_adjProlif, lab, colour = sc)) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.5), size = 0.2, linewidth = 0.4) + scale_x_log10() +
  scale_colour_manual(values = c(`7-gene score` = cN, `CCP score` = PAL[["accent"]]), name = NULL) + labs(x = "HR, log scale", y = NULL) + theme_pub() + theme(legend.position = "bottom", legend.margin = margin(-2, 0, 0, 0, "mm"))
fin(A13 | B13 | C13, "SuppFig_S13.pdf", 178, 62)

sc1 <- fread(file.path(P, "figures/celltype_localization/celltype_candidate_expression.csv"))
lin <- c(Luminal = "Epithelial", Basal = "Epithelial", `Proximal Luminal` = "Epithelial", Endothelial = "Endothelial", Fibroblast = "Stromal", `Smooth muscle cell` = "Stromal",
         Pericyte = "Stromal", SKMC = "Stromal", Schwann = "Stromal", CD4 = "Immune", CD8 = "Immune", B = "Immune", Plasma = "Immune", `NK/ILC` = "Immune", Neutrophil = "Immune",
         `Mast/Basophil` = "Immune", Macrophages_M0 = "Immune", Macrophages_M1 = "Immune", Macrophages_M2 = "Immune", Macrophage_Markers = "Immune")
sc1[, lineage := lin[celltype]]; stopifnot(!anyNA(sc1$lineage))
sc1[, z := as.numeric(scale(mean_expr)), by = .(dataset, gene)][is.na(z), z := 0]
GEN <- c("KIF18B", "FOXM1", "NEK2", "TOP2A")
ag <- sc1[gene %in% GEN, .(z = weighted.mean(z, n), pct = weighted.mean(pct_expr, n)), by = .(dataset, gene, lineage)]
ag[, `:=`(gene = factor(gene, levels = rev(GEN)), lineage = factor(lineage, levels = c("Epithelial", "Endothelial", "Stromal", "Immune")))]
note("S14: %d datasets, %s cells; genes %s; epithelial mean scaled expression by gene: %s", uniqueN(sc1$dataset), format(sc1[, sum(n), by = .(dataset, celltype)][, sum(V1)] / uniqueN(sc1$gene), big.mark = ","),
     paste(GEN, collapse = ", "), ag[lineage == "Epithelial", .(m = mean(z)), by = gene][, paste(sprintf("%s %.2f", gene, m), collapse = "; ")])
S14 <- ggplot(ag, aes(dataset, gene)) + geom_point(aes(size = pct, fill = z), shape = 21, colour = "grey35", stroke = 0.25) + facet_grid(~lineage) +
  scale_fill_gradient2(low = cN, mid = "white", high = cT, midpoint = 0, name = "Scaled\nexpression") + scale_size_area(max_size = 3.5, name = "Cells\nexpressing, %") +
  labs(x = NULL, y = NULL) + theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face = "italic"), plot.margin = margin(6, 2, 2, 2, "mm"))
save_panel(S14, file.path(O, "SuppFig_S14.pdf"), 178, 62)

ci <- fread(file.path(P, "nested_risk_model_results/Cindex.csv")); dc <- fread(file.path(P, "nested_risk_model_results/deltaC_bootstrap.csv")); lrt <- fread(file.path(P, "nested_risk_model_results/LR_tests.csv"))
ci[, lab := factor(c("Clinical", "+ Expression", "+ ASE")[seq_len(.N)], levels = c("Clinical", "+ Expression", "+ ASE"))]
note("S16A nested Cox (n=101, 28 events): corrected C-index %s; apparent %s; LR P %s; dC %s", ci[, paste(sprintf("%s %.3f", lab, c_corrected), collapse = ", ")],
     ci[, paste(sprintf("%.3f", c_apparent), collapse = ", ")], lrt[, paste(sprintf("%.2f", p), collapse = ", ")], dc[, paste(sprintf("%s %.3f [%.3f, %.3f]", comparison, dC_mean, ci_lo, ci_hi), collapse = "; ")])
A16 <- ggplot(ci, aes(lab, c_corrected)) + geom_point(size = 1.8) + geom_point(aes(y = c_apparent), shape = 21, size = 1.8, fill = "white") +
  geom_hline(yintercept = ci$c_corrected[1], linetype = 2, linewidth = 0.3) + scale_y_continuous(limits = c(0.5, 0.8)) + labs(x = NULL, y = "C-index") + theme_pub()
B16 <- ggplot(dc, aes(dC_mean, comparison)) + geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) + geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), size = 0.25, linewidth = 0.4) +
  scale_y_discrete(labels = c(`M2-M1` = "+ Expression", `M3-M2` = "+ ASE", `M3-M1` = "+ Expression + ASE")) + labs(x = "Change in C-index", y = NULL) + theme_pub()
cli <- fread("/path/to/data_root/CPGEA/CPGEA_web/CPGEA_Clinical_Info.tsv")[, .(pid = as.character(`Sequencing ID`), gs = suppressWarnings(as.numeric(GS)),
  bcr = fifelse(BCR == "YES", 1L, fifelse(BCR == "NO", 0L, NA_integer_)), months = suppressWarnings(as.numeric(sub("m$", "", `Month to surgery/BCR`))))]
pur <- fread(file.path(V, "delta_eqtl_clean/delta_clinical_patient_table.csv"))[, .(pid = as.character(pid), purity)]
cl <- merge(merge(fread(file.path(V, "delta_ase_clean/dilution_lowbaseline_patients.csv"))[, .(pid = as.character(pid), gain)], cli, by = "pid"), pur, by = "pid")
g <- cl[is.finite(gs)]; st <- cor.test(g$gain, g$gs, method = "spearman", exact = FALSE)
note("S16C: weak-baseline gain vs Gleason, n=%d: Spearman rho %.3f, P=%.2f", nrow(g), st$estimate, st$p.value)
note("S16C axis -0.02 to 0.04: patients outside the axis %d (%s)", g[gain < -0.02 | gain > 0.04, .N], g[gain < -0.02 | gain > 0.04, paste(sprintf("%.3f", gain), collapse = ", ")])
C16 <- ggplot(g, aes(factor(gs), gain)) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) + geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3, fill = "grey85") +
  geom_jitter(width = 0.15, height = 0, size = 0.4, alpha = 0.6) + coord_cartesian(ylim = c(-0.02, 0.04)) + labs(x = "Gleason score", y = "Weak-baseline gain") + theme_pub()
int <- function(x) qnorm((rank(x) - 0.5) / length(x))
dd <- cl[is.finite(bcr) & is.finite(months) & months > 0][, z := int(gain)]
HR <- rbindlist(lapply(list(Unadjusted = Surv(months, bcr) ~ z, `Adjusted, Gleason + purity` = Surv(months, bcr) ~ z + gs + purity), function(f) {
  m <- summary(coxph(f, dd)); data.table(hr = m$conf.int[1, 1], lo = m$conf.int[1, 3], hi = m$conf.int[1, 4], p = m$coef[1, 5]) }), idcol = "lab")
note("S16D BCR, weak-baseline gain, HR per SD of rank-normalized gain (n=%d, %d events): %s", nrow(dd), sum(dd$bcr), HR[, paste(sprintf("%s %.2f [%.2f, %.2f], P=%.2f", lab, hr, lo, hi, p), collapse = "; ")])
D16 <- ggplot(HR, aes(hr, factor(lab, levels = rev(lab)))) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) + geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.25, linewidth = 0.4) +
  labs(x = "BCR hazard ratio per SD", y = NULL) + theme_pub()
fin((A16 | B16) / (C16 | D16), "SuppFig_S16.pdf", 150, 100)
writeLines(LS, file.path(O, "supp_C_legend_stats.txt")); cat("DONE\n")
