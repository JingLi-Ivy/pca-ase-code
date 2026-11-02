#!/usr/bin/env Rscript

library(facets)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_facets.R <pileup_file> <output_prefix>")
}

pileup_file <- args[1]
output_prefix <- args[2]

cat("Reading pileup file:", pileup_file, "\n")
rcmat <- readSnpMatrix(pileup_file)

cat("Preprocessing data...\n")
xx <- preProcSample(rcmat)

cat("Running the facets algorithm...\n")
rets <- list()
dipLogR <- NULL

rets[[1]] <- procSample(xx, cval=150)
dipLogR <- rets[[1]]$dipLogR

rets[[2]] <- procSample(xx, cval=100, dipLogR=dipLogR)

fit <- procSample(xx, cval=50, dipLogR=dipLogR)

cat("Running EM estimation...\n")
emfit <- emcncf(fit)

cat("Saving results to:", paste0(output_prefix, "_facets_fits.rds"), "\n")
saveRDS(rets, file=paste0(output_prefix, "_facets_fits.rds"))
saveRDS(fit, file=paste0(output_prefix, "_facets_fit.rds"))
saveRDS(emfit, file=paste0(output_prefix, "_facets_emfit.rds"))

cat("Writing result summary...\n")
purity <- emfit$purity
ploidy <- emfit$ploidy
dipLogR <- fit$dipLogR
alBalLogR <- fit$alBalLogR
flags <- paste(fit$flags, collapse=", ")

cat("Saving summary to:", paste0(output_prefix, "_summary.txt"), "\n")
sink(paste0(output_prefix, "_summary.txt"))
cat("Sample:", basename(output_prefix), "\n")
cat("Purity:", purity, "\n")
cat("Ploidy:", ploidy, "\n")
cat("dipLogR:", dipLogR, "\n")
cat("alBalLogR:", alBalLogR, "\n")
cat("Flags:", flags, "\n")
sink()

cat("Summarizing CNV results...\n")
cncf <- emfit$cncf

cat("Saving cncf data to:", paste0(output_prefix, "_cncf.csv"), "\n")
write.csv(cncf, file=paste0(output_prefix, "_cncf.csv"), row.names=FALSE)

cnv_data <- cncf

high_cnv_segs <- cnv_data[cnv_data$tcn.em > 20, ]
num_high_cnv_segs <- nrow(high_cnv_segs)

if (num_high_cnv_segs > 0) {
  cat("Found", num_high_cnv_segs, "high-copy-number segments (CNV > 20)\n")
  write.table(high_cnv_segs, 
              file=paste0(output_prefix, "_high_cnv_segments.txt"), 
              sep="\t", 
              row.names=FALSE, 
              quote=FALSE)
} else {
  cat("No high-copy-number segments (CNV > 20) found\n")
}

cnv_counts <- table(cut(cnv_data$tcn.em, 
                        breaks=c(0, 1, 2, 3, 5, 10, 20, Inf), 
                        labels=c("0-1", "1-2", "2-3", "3-5", "5-10", "10-20", ">20"),
                        right=TRUE))

sink(paste0(output_prefix, "_cnv_statistics.txt"))
cat("CNV level distribution:\n")
print(cnv_counts)
cat("\nTotal segments:", nrow(cnv_data), "\n")
cat("High-copy-number segments (CNV > 20):", num_high_cnv_segs, "\n")
if (num_high_cnv_segs > 0) {
  cat("\nHigh-copy-number segment details:\n")
  print(high_cnv_segs[, c("chrom", "start", "end", "tcn.em", "lcn.em")])
}
sink()

cat("Drawing plots...\n")
pdf(paste0(output_prefix, "_facets_plots.pdf"), width=14, height=10)

plotSample(x=fit, emfit=emfit)

plotCNCF(fit, emfit=emfit)

dev.off()

cat("FACETS analysis done\n")
