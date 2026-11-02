#!/usr/bin/env Rscript
# 59a (2026-10-12): Fig 2A alone for v10. Panel code = lines 16-40 of 59_fig2_pub.R (ASE cohort only, 118 tumors), saved with its tag at the
# size of panel A in the v9-assembled Fig 2 (B-G unchanged), so that it can replace the 194-tumor panel drawn by 24_fig2A_fig7B.R.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
V <- "/path/to/data_root/ASE/code/v8_round2"; source(file.path(V, "fig_theme.R"))
P <- "/path/to/data_root/ASE"; O <- file.path(V, "pub_figs_v9"); dir.create(O, showWarnings = FALSE)
RES <- "/path/to/data_root/variants/batch_43/facets/results"
LS <- character(0); note <- function(...) { LS <<- c(LS, sprintf(...)); cat(sprintf(...), "\n") }
tags <- theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))

coh <- sub("_ase_results\\.csv$", "", list.files(P, pattern = "^T\\d+_ase_results\\.csv$"))
chrlen <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,
            133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)
off <- c(0, cumsum(as.numeric(chrlen)))[1:22]; names(off) <- 1:22; BIN <- 5e6
fs <- file.path(RES, paste0(coh, "_cncf.csv")); fs <- fs[file.exists(fs)]
bm <- rbindlist(lapply(fs, function(f) { d <- fread(f)[chrom %in% 1:22 & is.finite(cnlr.median) & is.finite(start) & is.finite(end)]
  if (!nrow(d)) return(NULL)
  d[, `:=`(b0 = floor((off[as.character(chrom)] + start) / BIN), b1 = floor((off[as.character(chrom)] + end) / BIN))]
  rbindlist(lapply(seq_len(nrow(d)), function(i) data.table(sample = sub("_cncf.csv", "", basename(f)), bin = d$b0[i]:d$b1[i], cnlr = d$cnlr.median[i]))) }))
bm <- bm[, .(cnlr = median(cnlr)), by = .(sample, bin)]
fga <- fread(file.path(P, "figures/supp_panels/facets_purity_ploidy_fga_per_sample.csv"))[sample %in% coh & is.finite(fga)][order(fga)]
bm <- bm[sample %in% fga$sample][, sample := factor(sample, levels = fga$sample)][, cnlr := pmax(pmin(cnlr, 1), -1)]
note("Fig2A: %d ASE-cohort tumors with FACETS segments (of %d ASE tumors); %d bins of 5 Mb; median FGA %.3f", uniqueN(bm$sample), length(coh), uniqueN(bm$bin), median(fga$fga))
cb <- floor(off / BIN)
A2 <- ggplot(bm, aes(bin, sample, fill = cnlr)) + geom_tile(width = 1, height = 1) +
  scale_fill_gradient2(low = PAL[["normal"]], mid = "white", high = PAL[["tumor"]], midpoint = 0, limits = c(-1, 1), name = "Log ratio") +
  geom_vline(xintercept = cb[-1] - 0.5, colour = "grey60", linewidth = 0.15) +
  scale_x_continuous(breaks = (cb + c(cb[-1], max(bm$bin)))/2, labels = ifelse(1:22 %% 2 == 1, 1:22, ""), expand = c(0, 0)) +
  labs(x = "Chromosome", y = sprintf("Tumors (n = %d)", uniqueN(bm$sample))) + theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), legend.key.height = unit(4, "mm"))

save_panel(A2 + labs(tag = "A") + tags, file.path(O, "Fig2A_118.pdf"), 129, 71)
writeLines(LS, file.path(O, "fig2A_118_legend_stats.txt"))
