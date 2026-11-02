#!/usr/bin/env Rscript
# 59b (2026-10-13): publication redraw of main Fig 2 for v10, in the shared fig_theme.R style.
# Replaces the v9 collage, whose panels B-G came from older scripts in mixed styles (overlapping and squeezed text).
# Panel content follows the v10 legend; the pie chart (old D) is dropped because the UpSet plot carries the same
# exclusive combinations, so old E/F/G become D/E/F:
#  A copy-number landscape of the 118 ASE-cohort tumors (built by 59a_fig2A_panel.R, reused unchanged)
#  B SCNA-ASE Spearman correlation per PCa-ASE gene along the genome (Supplemental Table S3)
#  C deletion- versus amplification-associated genes (FDR < 0.1)
#  D gene-level source attribution as an UpSet plot: bars = mutually exclusive evidence combinations,
#    left bars = non-exclusive totals (Supplemental Table S5); unassigned genes stated in the panel
#  E, F ATP2C1 and MCM9 examples (per-tumor max |a| versus SCNA score, as 04_Fig2_ASE_CNV_analysis.R)
# Copy-number colours: loss blue, gain red (as A). ASE status is shown in greys so that red/blue keep one meaning.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
P <- "/path/to/data_root/ASE"; C <- file.path(P, "code"); O <- file.path(V, "pub_figs_v10"); dir.create(O, showWarnings = FALSE)
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tg <- function(p, t) p + labs(tag = t) + theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
sci <- function(x, d = 1) { e <- floor(log10(x)); m <- x / 10^e
  if (round(m, d) == 10) { m <- 1; e <- e + 1 }; sprintf("%s %%*%% 10^%d", formatC(m, format = "f", digits = d), e) }
LOSS <- PAL[["normal"]]; GAIN <- PAL[["tumor"]]

ea <- new.env(); sys.source(file.path(V, "59a_fig2A_panel.R"), envir = ea); A2 <- ea$A2

cc <- fread(file.path(C, "fig2_results/pca_ase_cnv_correlations.csv"))[is.finite(rho) & is.finite(pval)]
gtf <- fread(cmd = "zcat /path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz | awk -F'\\t' '$3==\"gene\"' | cut -f1,4,9", header = FALSE, sep = "\t", quote = "")
gtf[, gene_id_clean := sub('.*gene_id "([^"]+)".*', "\\1", V3)]
cc <- merge(cc, gtf[V1 %in% 1:22, .(gene_id_clean, chrom = as.integer(V1), pos = V2)], by = "gene_id_clean")
chrlen <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,
            133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)
off <- c(0, cumsum(as.numeric(chrlen)))[1:22]
cc[, x := off[chrom] + pos]
cc[, cl := factor(fifelse(fdr < 0.1, fifelse(rho < 0, "Deletion-associated", "Amplification-associated"), "Not significant"),
                  levels = c("Deletion-associated", "Amplification-associated", "Not significant"))]
pthr <- cc[fdr < 0.1, max(pval)]
note("Fig2B: %d PCa-ASE genes with SCNA-ASE correlation (mapped to autosomes %d); FDR<0.1: %d (deletion %d, amplification %d); P threshold %.3g",
     nrow(fread(file.path(C, "fig2_results/pca_ase_cnv_correlations.csv"))), nrow(cc), cc[fdr < 0.1, .N], cc[cl == "Deletion-associated", .N], cc[cl == "Amplification-associated", .N], pthr)
colsC <- c(`Deletion-associated` = LOSS, `Amplification-associated` = GAIN, `Not significant` = "grey75")
mid <- off + chrlen / 2
B2 <- ggplot(cc[order(cl == "Not significant", decreasing = TRUE)], aes(x, -log10(pval), colour = cl)) +
  geom_vline(xintercept = off[-1], colour = "grey90", linewidth = 0.2) +
  geom_point(size = 0.5) + geom_hline(yintercept = -log10(pthr), linetype = 2, linewidth = 0.3) +
  scale_colour_manual(values = colsC, name = NULL) +
  scale_x_continuous(breaks = mid, labels = ifelse(1:22 %% 2 == 1, 1:22, ""), expand = c(0.005, 0)) +
  labs(x = "Chromosome", y = expression(-log[10]~italic(P))) + theme_pub() +
  theme(axis.ticks.x = element_blank(), legend.position = "top", legend.justification = "right",
        legend.margin = margin(0, 0, -2, 0, "mm"), legend.key.size = unit(2.5, "mm")) + guides(colour = guide_legend(override.aes = list(size = 1.5)))
cd <- cc[cl != "Not significant", .N, by = cl][, `:=`(cl = droplevels(cl), pct = 100 * N / sum(N))]
C2 <- ggplot(cd, aes(cl, N, fill = cl)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", N, pct)), vjust = -0.4, size = BASE / .pt, family = FONT) +
  scale_fill_manual(values = colsC, guide = "none") + scale_y_continuous(limits = c(0, 700), expand = c(0, 0)) +
  scale_x_discrete(labels = c(`Deletion-associated` = "Deletion-\nassociated", `Amplification-associated` = "Amplification-\nassociated")) +
  labs(x = NULL, y = "SCNA-associated genes") + theme_pub()

g5 <- fread(file.path(C, "mechanism_integration_v2/full_integration_results.csv"))
S <- c(SCNA = "is_SCNA", eQTL = "is_eQTL", sQTL = "is_sQTL", Methylation = "is_Methylation")
M <- as.data.table(lapply(S, function(s) g5[[s]] %in% TRUE))
M[, combo := apply(.SD, 1, function(r) paste(names(S)[r], collapse = "+"))]
cmb <- M[combo != "", .N, by = combo][order(-N)][, id := factor(combo, levels = combo)]
nun <- M[combo == "", .N]; setsz <- data.table(set = factor(names(S), levels = rev(names(S))), n = sapply(names(S), function(s) sum(M[[s]])))
note("Fig2D: %d PCa-ASE genes; non-exclusive %s; unassigned %d (%.1f%%); exclusive combinations: %s", nrow(g5),
     setsz[, paste(sprintf("%s %d (%.1f%%)", set, n, 100 * n / nrow(g5)), collapse = ", ")], nun, 100 * nun / nrow(g5), cmb[, paste(sprintf("%s %d", combo, N), collapse = "; ")])
mat <- CJ(id = cmb$id, set = names(S))[, on := mapply(function(i, s) s %in% strsplit(as.character(i), "+", fixed = TRUE)[[1]], id, set)]
mat[, set := factor(set, levels = rev(names(S)))]
Dbar <- ggplot(cmb, aes(id, N)) + geom_col(fill = "grey35", width = 0.65) +
  geom_text(aes(label = N), vjust = -0.4, size = (BASE - 1) / .pt, family = FONT) +
  annotate("text", x = nrow(cmb), y = 470, hjust = 1, vjust = 1, family = FONT, size = BASE / .pt,
           label = sprintf("Unassigned: %s genes (%.1f%%)", format(nun, big.mark = ","), 100 * nun / nrow(g5))) +
  scale_y_continuous(limits = c(0, 520), expand = c(0, 0)) + labs(x = NULL, y = "Genes (exclusive)") + theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
seg <- mat[on == TRUE, .(lo = min(as.integer(set)), hi = max(as.integer(set))), by = id]
Dmat <- ggplot(mat, aes(id, set)) + geom_point(aes(colour = on), size = 1.8) +
  geom_segment(data = seg, aes(x = id, xend = id, y = lo, yend = hi), inherit.aes = FALSE, linewidth = 0.5, colour = "grey20") +
  scale_colour_manual(values = c(`TRUE` = "grey20", `FALSE` = "grey88"), guide = "none") + labs(x = NULL, y = NULL) + theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
Dset <- ggplot(setsz, aes(n, set)) + geom_col(fill = "grey55", width = 0.6) + scale_x_reverse(breaks = c(600, 300, 0), expand = c(0, 0)) +
  labs(x = "Genes (total)", y = NULL) + theme_pub() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
D2 <- (plot_spacer() + Dbar + Dset + Dmat) + plot_layout(ncol = 2, widths = c(0.28, 1), heights = c(1, 0.55))

ex <- c(ATP2C1 = "ENSG00000017260", MCM9 = "ENSG00000111877")
af <- list.files(P, "^T\\d+_ase_results\\.csv$", full.names = TRUE)
ase <- rbindlist(lapply(af, function(f) { d <- fread(f)[sub("\\..*", "", gene_id) %chin% ex]; if (!nrow(d)) return(NULL)
  d[, .(Sample = sub("_ase_results.csv", "", basename(f)), gid = sub("\\..*", "", gene_id), a = abs(a), fdr)] }))
ase <- ase[, .(a = max(a, na.rm = TRUE), sig = any(fdr < 0.1, na.rm = TRUE)), by = .(Sample, gid)][is.finite(a)]
cn <- fread(cmd = "awk -F'\\t' 'NR==1 || $9==\"ATP2C1\" || $9==\"MCM9\"' /path/to/data_root/variants/batch_43/facets/annotated_cnv_output.tsv")
cn <- cn[, .(cnv = mean(log2(pmax(tcn, 0.1) / 2))), by = .(Sample = sample_id, gene = gene_name)][, gid := ex[gene]]
ex_plot <- function(g) {
  d <- merge(ase[gid == ex[[g]]], cn[gene == g], by = c("Sample", "gid"))
  r <- cc[gene_id_clean == ex[[g]]]; ct <- cor.test(d$a, d$cnv, method = "spearman", exact = FALSE)
  note("Fig2 %s: n=%d tumors; recomputed rho %.3f (S3 rho %.3f, FDR %.2g); ASE-significant tumors %d", g, nrow(d), ct$estimate, r$rho, r$fdr, d[sig == TRUE, .N])
  lab <- sprintf("rho == %.2f*','~~FDR == %s", r$rho, sci(r$fdr))
  ggplot(d[order(sig)], aes(a, cnv, colour = sig)) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_point(size = 1.1, alpha = 0.85) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"), labels = c(`TRUE` = "ASE (FDR < 0.1)", `FALSE` = "No ASE"), name = NULL) +
    labs(x = expression("ASE magnitude, |"*italic(a)*"|"), y = "SCNA score", title = g, subtitle = parse(text = lab)[[1]]) +
    theme_pub() + theme(plot.title = element_text(size = BASE, face = "italic", family = FONT, margin = margin(0, 0, 0.5, 0, "mm")),
                        plot.subtitle = element_text(size = BASE, family = FONT, margin = margin(0, 0, 1, 0, "mm"))) }
E2 <- ex_plot("ATP2C1"); F2 <- ex_plot("MCM9")
EF <- (tg(E2, "E") / tg(F2, "F")) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

g <- wrap_elements(full = tg(A2, "A")) / (tg(B2, "B") + tg(C2, "C") + plot_layout(widths = c(2.1, 1))) /
     (wrap_elements(full = D2) + labs(tag = "D") + theme(plot.tag = element_text(size = 10, face = "bold", family = FONT)) + EF + plot_layout(widths = c(2.1, 1))) +
     plot_layout(heights = c(1.05, 0.95, 1.35))
save_panel(g, file.path(O, "Fig2.pdf"), 178, 205)
writeLines(LS, file.path(O, "fig2_legend_stats.txt")); cat("DONE\n")
