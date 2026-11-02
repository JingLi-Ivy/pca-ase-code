#!/usr/bin/env Rscript
# A13: KIF18B / KIF18B-DT bidirectional-promoter locus figure (hg38).
# A overview: gene models, ASE-reporting SNPs, LD (r2) with rs9914422 in CPGEA, ENCODE cCREs.
# B promoter zoom: TSSs, CpG island, cCREs, reporter fragment, ENCODE TF ChIP-seq clusters (selected factors).
# C base level: sequence +/-10 bp, both alleles, best-scoring JASPAR sites (NRF1 on C; SP2/KLF5/ZNF682/YY1 on G).
# D motif score by allele (JASPAR 2024).  E ADASTRA v6.1 allele-specific binding at rs9914422.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(jsonlite); library(ggrepel)})
D <- "/path/to/data_root/ASE/code/v8_round2/locus_fig"; setwd(D)
RS <- 44947872; FR <- c(44947189, 44948553); OV <- c(44918000, 44960000); ZM <- c(44946900, 44949000)
th <- theme_classic(base_size = 9) + theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1), legend.key.size = unit(0.3, "cm"))
xs <- function(lim) if (diff(lim) > 1e4) scale_x_continuous(limits = lim, expand = c(0, 0), labels = function(x) sprintf("%.2f", x / 1e6), name = "chr17 (Mb, hg38)") else
  scale_x_continuous(limits = lim, expand = c(0, 0), labels = function(x) format(x, big.mark = ","), name = "chr17 position (hg38)")

ex <- fread("gene_exons.tsv", header = FALSE, col.names = c("gene","tx","type","start","end","strand"))
ex <- ex[type == "exon"]
txs <- ex[, .(start = min(start), end = max(end), strand = strand[1], gene = gene[1]), by = tx]
txs <- txs[order(gene, tx)][, y := -.I]
ex <- merge(ex, txs[, .(tx, y)], by = "tx")
pg <- function(lim) ggplot() + geom_segment(data = txs, aes(x = pmax(start, lim[1]), xend = pmin(end, lim[2]), y = y, yend = y), linewidth = 0.3) +
  geom_rect(data = ex, aes(xmin = pmax(start, lim[1]), xmax = pmin(end, lim[2]), ymin = y - 0.3, ymax = y + 0.3, fill = gene)) +
  geom_text(data = txs[, .(x = ifelse(gene == "KIF18B-DT", pmin(end, lim[2]), pmax(start, lim[1])), h = as.numeric(gene == "KIF18B-DT"), y, lab = paste0(gene, " ", tx, ifelse(strand == "-", " (-)", " (+)")))], aes(x, y + 0.5, label = lab, hjust = h), size = 2) +
  scale_fill_manual(values = c(KIF18B = "#2166AC", `KIF18B-DT` = "#B2182B"), guide = "none") + geom_vline(xintercept = RS, colour = "red", linetype = 2, linewidth = 0.3) +
  xs(lim) + labs(y = "Genes") + th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ld <- fread("ld_rs9914422.vcor"); ld <- ld[, .(pos = POS_B, r2 = UNPHASED_R2, id = ID_B)]; ld <- rbind(ld, data.table(pos = RS, r2 = 1, id = "rs9914422"))
ase <- fread("../kif18b_audit/kif18b_rna_snps.csv")[gene_id == "ENSG00000186185", .(tumors = uniqueN(sample_id), reads = sum(N)), by = .(pos = start)]
pl <- ggplot(ld[pos >= OV[1] & pos <= OV[2]], aes(pos, r2, colour = r2)) + geom_point(size = 0.9) +
  geom_point(data = ld[id == "rs9914422"], colour = "red", size = 2.2, shape = 18) +
  geom_text_repel(data = ld[id == "rs9914422"], aes(label = id), size = 2.3, colour = "red", nudge_y = -0.25, min.segment.length = 0) +
  geom_rug(data = ase[pos >= OV[1] & pos <= OV[2]], aes(x = pos), inherit.aes = FALSE, sides = "b", colour = "darkgreen", length = unit(0.06, "npc")) +
  scale_colour_gradient(low = "grey75", high = "navy", limits = c(0, 1), name = expression(r^2)) + xs(OV) +
  labs(y = expression(r^2 ~ "with rs9914422")) + th
cc <- as.data.table(fromJSON("ccre.json")); lab <- if ("ucscLabel" %in% names(cc)) "ucscLabel" else "encodeLabel"
cc[, cls := get(lab)]
pc <- function(lim) ggplot(cc[chromEnd >= lim[1] & chromStart <= lim[2]]) + geom_rect(aes(xmin = pmax(chromStart, lim[1]), xmax = pmin(chromEnd, lim[2]), ymin = 0, ymax = 1, fill = cls)) +
  scale_fill_manual(values = c(prom = "#FF0000", PLS = "#FF0000", enhP = "#FFA700", pELS = "#FFA700", enhD = "#FFCD00", dELS = "#FFCD00", K4m3 = "#FFAAAA", CTCF = "#00B0F0"), name = "ENCODE cCRE") +
  geom_vline(xintercept = RS, colour = "red", linetype = 2, linewidth = 0.3) + xs(lim) + labs(y = "cCRE") + th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
eq <- fromJSON("gtex_eqtl.json")$rs144461621
nld <- ld[r2 >= 0.8 & id != "rs9914422" & pos >= OV[1] & pos <= OV[2], .N]
A <- (pg(OV) / pl / pc(OV)) + plot_layout(heights = c(1.4, 1.2, 0.35)) +
  plot_annotation(title = sprintf("A  KIF18B locus: %d other CPGEA SNPs in r2 >= 0.8 with rs9914422; green ticks = %d ASE-reporting SNPs.\n    GTEx prostate lead KIF18B eQTL rs144461621 lies %.0f kb away (chr17:%s), outside this window.", nld, nrow(ase[pos >= OV[1] & pos <= OV[2]]), (RS - eq) / 1e3, format(eq, big.mark = ",")),
                  theme = theme(plot.title = element_text(size = 9)))

tf <- as.data.table(fromJSON("encode_tfbs_clusters.json")$encRegTfbsClustered)
SEL <- c("FOXM1","MYBL2","E2F4","E2F1","GABPA","ELF1","ETV6","NRF1","SP1","YY1","FOXK2")
tfz <- tf[name %in% SEL & chromEnd >= ZM[1] & chromStart <= ZM[2]][, name := factor(name, levels = rev(SEL))]
tss <- data.table(x = c(44947773, 44947912), lab = c("KIF18B TSS (-)", "KIF18B-DT TSS (+)"), dir = c(-1, 1))
cpg <- data.table(s = 44947290, e = 44947976)
B1 <- ggplot() + annotate("rect", xmin = cpg$s, xmax = cpg$e, ymin = 2.2, ymax = 2.8, fill = "darkolivegreen3") + annotate("text", x = cpg$s - 20, y = 2.5, label = "CpG island", size = 2.2, hjust = 1) +
  annotate("rect", xmin = FR[1], xmax = FR[2], ymin = 1.2, ymax = 1.8, fill = "grey70") + annotate("text", x = FR[2] - 20, y = 1.95, label = "luciferase reporter fragment (1,364 bp)", size = 2.2, hjust = 1) +
  geom_segment(data = tss, aes(x = x, xend = x + dir * 180, y = 0.5, yend = 0.5), arrow = arrow(length = unit(0.12, "cm")), colour = c("#2166AC", "#B2182B")) +
  geom_text(data = tss, aes(x = x + dir * 200, y = 0.5, label = lab, hjust = ifelse(dir > 0, 0, 1)), size = 2.1) +
  geom_vline(xintercept = RS, colour = "red", linetype = 2, linewidth = 0.3) + annotate("text", x = RS + 20, y = 3.3, label = "rs9914422 (G>C)", colour = "red", size = 2.3, hjust = 0) +
  xs(ZM) + ylim(0, 3.6) + labs(y = NULL) + th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
B2 <- ggplot(tfz) + geom_rect(aes(xmin = pmax(chromStart, ZM[1]), xmax = pmin(chromEnd, ZM[2]), ymin = as.numeric(name) - 0.35, ymax = as.numeric(name) + 0.35, alpha = score), fill = "black") +
  scale_y_continuous(breaks = seq_along(levels(tfz$name)), labels = levels(tfz$name)) + scale_alpha(range = c(0.25, 1), guide = "none") +
  geom_vline(xintercept = RS, colour = "red", linetype = 2, linewidth = 0.3) + xs(ZM) + labs(y = "ENCODE\nTF ChIP") + th
B <- (B1 / pc(ZM) / B2) + plot_layout(heights = c(1, 0.3, 1.6)) + plot_annotation(title = "B  Bidirectional promoter (ENCODE TF ChIP-seq clusters; 150 factors overlap rs9914422)", theme = theme(plot.title = element_text(size = 9)))

seq <- toupper(gsub("\\s", "", paste(readLines("seq_pm30.txt", warn = FALSE), collapse = "")))
off <- -12:12; bases <- strsplit(substr(seq, 31 - 12, 31 + 12), "")[[1]]
ms <- fread("motif_best_sites.tsv")
show <- rbind(ms[tf == "NRF1" & allele == "C"], ms[tf %in% c("SP2","KLF5","ZNF682","YY1") & allele == "G"])
show[, `:=`(row = .I, a = pmax(start1 - RS, -12), b = pmin(end1 - RS, 12))]
bt <- data.table(x = off, G = bases, C = replace(bases, 13, "C"))
C1 <- ggplot() + geom_text(data = bt, aes(x, 2, label = G, colour = x == 0), family = "mono", size = 3.2) +
  geom_text(data = bt, aes(x, 1.2, label = C, colour = x == 0), family = "mono", size = 3.2) +
  annotate("text", x = -13, y = c(2, 1.2), label = c("G allele", "C allele"), hjust = 1, size = 2.4) +
  geom_rect(data = show, aes(xmin = a - 0.45, xmax = b + 0.45, ymin = 0.35 - row * 0.55, ymax = 0.7 - row * 0.55, fill = allele)) +
  geom_text(data = show, aes(x = 13.2, y = 0.525 - row * 0.55, label = sprintf("%s %s (%s strand), best on %s allele: %.2f", tf, jaspar, strand, allele, rel_score)), hjust = 0, size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "black"), guide = "none") +
  scale_fill_manual(values = c(C = "#B2182B", G = "grey40"), guide = "none") +
  scale_x_continuous(limits = c(-19, 32), breaks = c(-10, -5, 0, 5, 10), name = "bp relative to rs9914422 (+ strand, hg38)") +
  labs(y = NULL, title = "C  Sequence context and best JASPAR 2024 sites (red box: site strengthened by C; grey: site weakened by C)") + th +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), plot.title = element_text(size = 9))

sc <- as.data.table(fromJSON("motif_allele_scores.json")); if (ncol(sc) == 5) setnames(sc, c("tf","id","G","C","d")) else sc <- sc[, 1:5]
sc <- sc[!is.na(suppressWarnings(as.numeric(G)))][, `:=`(G = as.numeric(G), C = as.numeric(C), d = as.numeric(d))]
sc[, tf := factor(tf, levels = sc[order(d), tf])]
Dp <- ggplot(sc) + geom_segment(aes(x = G, xend = C, y = tf, yend = tf), colour = "grey60") + geom_point(aes(G, tf), colour = "grey30", size = 1.6) +
  geom_point(aes(C, tf), colour = "#B2182B", size = 1.6) + geom_vline(xintercept = 0.8, linetype = 3) +
  labs(x = "relative motif score (grey G, red C)", y = NULL, title = "D  Motif score by allele") + th + theme(plot.title = element_text(size = 9))

ad <- fread("adastra_rs9914422.tsv")
ad <- ad[order(minfdr)][minfdr < 0.1 | tf %in% c("FOXM1","FOXA1","AR","MYBL2","E2F4")]
ad[, `:=`(lab = sprintf("%s  %s  FDR %.1e", tf, pref, minfdr), tf = factor(tf, levels = rev(tf)))]
Ep <- ggplot(ad, aes(y = tf)) + geom_col(aes(x = -log10(minfdr), fill = pref), width = 0.7) + geom_vline(xintercept = -log10(0.05), linetype = 3) +
  scale_fill_manual(values = c(G = "grey40", C = "#B2182B"), name = "preferred allele") +
  labs(x = expression(-log[10] ~ "FDR (allele-specific binding)"), y = NULL, title = "E  ADASTRA v6.1 at rs9914422") + th + theme(plot.title = element_text(size = 9))

pdf("KIF18B_locus_figure.pdf", width = 8.5, height = 14)
print(wrap_elements(A) / wrap_elements(B) / C1 / (Dp | Ep) + plot_layout(heights = c(1.25, 1.2, 0.5, 1)))
dev.off()
png("KIF18B_locus_figure.png", width = 8.5, height = 14, units = "in", res = 200)
print(wrap_elements(A) / wrap_elements(B) / C1 / (Dp | Ep) + plot_layout(heights = c(1.25, 1.2, 0.5, 1)))
dev.off()
cat("done\n")
