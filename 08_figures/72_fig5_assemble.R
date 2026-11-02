#!/usr/bin/env Rscript
# 72 (2026-10-12): main Fig 5 (v10), panels A-F assembled in R at 178 mm width (F = revised model, drawn here; replaces the Illustrator step).
# Panels are the ggplot objects of the scripts that own them (each script is run in its own environment; their own outputs are rewritten unchanged):
#   A = locus_panel() below (KIF18B-KIF18B-DT promoter, browser-track style; same inputs as 60 C5)
#   B = 25 A (phased, tumor-DNA-corrected C-haplotype fraction vs FOXM1; 0.486 [0.395, 0.579])
#   C = 60 F5 (dual-luciferase reporter)          D = 70 D5 (FOXM1 -> KIF18B, 118 tumors)          E = 60 G5 (siFOXM1 qRT-PCR)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
run <- function(f) { e <- new.env(); wd <- getwd(); on.exit(setwd(wd)); sys.source(file.path(V, f), envir = e); e }
e60 <- run("60_fig5_pub.R"); e25 <- run("25_kif18b_phase_panel.R"); e70 <- run("70_fig5d_suppS17.R")
locus_panel <- function() {
  L <- file.path(V, "locus_fig"); RS <- 44947872; FR <- c(44947189, 44948553); ZM <- c(44946900, 44949000); TS <- BASE / .pt
  gx <- fread(file.path(L, "gene_exons.tsv"), header = FALSE, col.names = c("gene", "tx", "type", "start", "end", "strand"))
  ge <- gx[type == "gene"]; ex <- unique(gx[type == "exon" & end >= ZM[1] & start <= ZM[2], .(gene, s = pmax(start, ZM[1]), e = pmin(end, ZM[2]))])
  ex <- ex[, .(s = min(s), e = max(e)), by = gene]
  cc <- as.data.table(jsonlite::fromJSON(file.path(L, "ccre.json")))[chromEnd >= ZM[1] & chromStart <= ZM[2]]
  cc[, `:=`(s = pmax(chromStart, ZM[1]), e = pmin(chromEnd, ZM[2]), cls = fifelse(encodeLabel == "PLS", "Promoter-like signature", "Proximal enhancer-like signature"))]
  yG <- c(KIF18B = 4.0, `KIF18B-DT` = 3.0); yR <- 1.9; yC <- 0.9; h <- 0.18
  gl <- ge[, .(gene, y = yG[gene], s = pmax(start, ZM[1]), e = pmin(end, ZM[2]), strand, tss = fifelse(strand == "-", end, start))]
  chev <- gl[, .(x = seq(s + 120, e - 120, by = 260)), by = .(gene, y, strand)][, `:=`(x0 = fifelse(strand == "-", x + 35, x - 35), x1 = x)]
  chev <- chev[!(gene == "KIF18B" & x > 44947600)]
  tss <- gl[, .(gene, y, tss, dir = fifelse(strand == "-", -1, 1))]
  note("Fig5A (redrawn): window chr17:%s-%s; exons %s; reporter %s-%s (%d bp); cCREs %d (%s)", format(ZM[1], big.mark = ","), format(ZM[2], big.mark = ","),
       ex[, paste(sprintf("%s %d-%d", gene, s, e), collapse = "; ")], format(FR[1], big.mark = ","), format(FR[2], big.mark = ","), diff(FR) + 1, nrow(cc), paste(cc$encodeLabel, collapse = ","))
  ggplot() +
    annotate("segment", x = RS, xend = RS, y = 0.5, yend = 4.75, colour = PAL[["tumor"]], linewidth = 0.35, linetype = "22") +
    annotate("text", x = RS, y = 4.85, label = "rs9914422", colour = PAL[["tumor"]], size = TS, family = FONT, vjust = 0) +
    geom_segment(data = gl, aes(x = s, xend = e, y = y, yend = y), linewidth = 0.35, colour = "grey30") +
    geom_segment(data = chev, aes(x = x0, xend = x1, y = y, yend = y), linewidth = 0.3, colour = "grey30", arrow = arrow(length = unit(0.9, "mm"), angle = 35)) +
    geom_rect(data = ex, aes(xmin = s, xmax = e, ymin = yG[gene] - h, ymax = yG[gene] + h), fill = "grey20") +
    geom_segment(data = tss, aes(x = tss, xend = tss, y = y + h, yend = y + 0.42), linewidth = 0.35) +
    geom_segment(data = tss, aes(x = tss, xend = tss + dir * 110, y = y + 0.42, yend = y + 0.42), linewidth = 0.35, arrow = arrow(length = unit(0.9, "mm"), angle = 35)) +
    annotate("text", x = c(44947200, 44948425), y = yG + 0.42, label = c("KIF18B", "KIF18B-DT"), fontface = "italic", size = TS, family = FONT, vjust = 0.5) +
    annotate("rect", xmin = FR[1], xmax = FR[2], ymin = yR - h, ymax = yR + h, fill = "grey60") +
    annotate("text", x = (FR[1] + RS) / 2, y = yR, label = sprintf("%s bp", format(diff(FR) + 1, big.mark = ",")), size = TS * 0.9, family = FONT, colour = "white") +
    geom_rect(data = cc, aes(xmin = s, xmax = e, ymin = yC - h, ymax = yC + h, fill = cls)) +
    scale_fill_manual(values = c(`Promoter-like signature` = "#E41A1C", `Proximal enhancer-like signature` = "#FFA700"), name = NULL) +
    scale_x_continuous(limits = ZM, expand = c(0, 0), breaks = seq(44947000, 44948500, 500), labels = function(x) format(x, big.mark = ",")) +
    scale_y_continuous(limits = c(0.45, 5.1), breaks = c(mean(yG), yR, yC), labels = c("Genes", "Reporter\nfragment", "ENCODE cCRE"), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + labs(x = "chr17 position (GRCh38)", y = NULL) + theme_pub() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 1, colour = "black"),
          legend.position = "bottom", legend.margin = margin(-1, 0, 0, 0, "mm"), legend.key.size = unit(2.5, "mm"))
}
model_panel <- function() {
  TS <- BASE / .pt; red <- PAL[["tumor"]]
  bx <- data.table(id = c("tp", "fx", "kf", "ph", "al"), x = c(18, 60, 105, 153, 105), y = c(40, 40, 40, 40, 13),
                   w = c(30, 20, 22, 40, 34), h = 9, fill = c("white", "#F4CACA", "#CFE3F3", "grey92", "white"), lt = c("22", "solid", "solid", "solid", "solid"),
                   lab = c("italic('TP53')*' / '*italic('RB1')*' loss'", "'FOXM1'", "italic('KIF18B')", "'Proliferation, migration'", "'rs9914422 C allele'"))
  ar <- data.table(x = c(33, 70, 116, 105), y = c(40, 40, 40, 17.5), xe = c(50, 94, 133, 105), ye = c(40, 40, 40, 35.5), lt = c("22", "solid", "solid", "solid"),
                   lab = c("published", "co-expression,\nsiRNA", "siRNA", "reporter: C > G"),
                   lx = c(41.5, 82, 124.5, 106.5), ly = c(42.8, 42.8, 42.8, 26.5), hj = c(0.5, 0.5, 0.5, 0))
  ggplot() +
    geom_rect(data = bx, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill, linetype = lt), colour = "grey25", linewidth = 0.35) +
    geom_text(data = bx, aes(x, y, label = lab), parse = TRUE, size = TS, family = FONT) +
    geom_segment(data = ar, aes(x, y, xend = xe, yend = ye, linetype = lt), linewidth = 0.45, arrow = arrow(length = unit(1.4, "mm"), type = "closed", angle = 25)) +
    geom_text(data = ar, aes(lx, ly, label = lab, hjust = hj), size = TS, family = FONT, vjust = 0, lineheight = 0.9) +
    annotate("curve", x = 60, y = 35.5, xend = 87.5, yend = 13, curvature = 0.25, colour = red, linetype = "11", linewidth = 0.55,
             arrow = arrow(length = unit(1.4, "mm"), type = "closed", angle = 25)) +
    annotate("text", x = 66, y = 21.5, label = "?", colour = red, size = TS * 1.6, family = FONT, fontface = "bold") +
    annotate("text", x = 58, y = 13, label = "C-selective amplification\nnot detected in tumors", colour = red, size = TS, family = FONT, hjust = 1, lineheight = 0.9) +
    annotate("segment", x = c(2, 62, 118), xend = c(12, 72, 128), y = 1.5, yend = 1.5, linetype = c("solid", "22", "11"), colour = c("black", "black", red), linewidth = 0.45) +
    annotate("text", x = c(14, 74, 130), y = 1.5, label = c("Supported", "Published", "Tested, not supported"), hjust = 0, size = TS, family = FONT) +
    scale_fill_identity() + scale_linetype_identity() + coord_cartesian(xlim = c(-7, 178), ylim = c(-1, 48), expand = FALSE, clip = "off") + theme_void()
}
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
pA <- locus_panel()
pB <- e25$A + labs(x = expression(italic("FOXM1")*" expression (FPKM)")); pC <- e60$F5; pD <- e70$D5; pE <- e60$G5
pF <- wrap_elements(full = model_panel())
g <- pA / ((free(pB, "label") | pC) + plot_layout(widths = c(1, 1))) / ((free(pD, "label") | pE) + plot_layout(widths = c(0.8, 1.2))) / pF +
  plot_layout(heights = c(48, 52, 50, 50)) + plot_annotation(tag_levels = "A") & tags
save_panel(g, file.path(V, "pub_figs_v9", "Fig5.pdf"), 178, 212)
writeLines(c(LS, readLines(file.path(V, "pub_figs_v9", "fig5_legend_stats.txt")), readLines(file.path(V, "kif18b_phase", "Fig5_KIF18B_phased_legend_stats.txt")),
             grep("^Fig5D", readLines(file.path(V, "pub_figs_v9", "fig5d_suppS17_legend_stats.txt")), value = TRUE)),
           file.path(V, "pub_figs_v9", "Fig5_legend_stats_sources.txt"))
cat("DONE\n")
