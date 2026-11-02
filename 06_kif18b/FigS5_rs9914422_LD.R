#!/usr/bin/env Rscript
# Supplementary Figure 5: rs9914422 KIF18B-locus haplotype LD structure.
# Pairwise LD (r^2) among the five FOXM1-modified haplotype SNPs + gene-body
# ASE-reporter rs8070376, computed from CPGEA germline genotypes (combined.134.qc).
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

PLINK <- "/path/to/data_root/CPGEA/CPGEA_genotype/combined.134.qc"
OUT   <- "/path/to/data_root/ASE/figures/rs9914422_LD"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

bim <- fread(paste0(PLINK, ".bim"))
fam <- fread(paste0(PLINK, ".fam"), header = FALSE)

read_snp <- function(rsid) {
  idx <- which(bim$V2 == rsid); if (length(idx) != 1) return(NULL)
  n <- nrow(fam); bps <- ceiling(n/4)
  con <- file(paste0(PLINK, ".bed"), "rb"); on.exit(close(con))
  readBin(con, "raw", 3); seek(con, 3 + (idx-1)*bps, origin = "start")
  bv <- as.integer(readBin(con, "raw", bps))
  bits <- as.integer(rbind(bv%%4, (bv%/%4)%%4, (bv%/%16)%%4, (bv%/%64)%%4))[seq_len(n)]
  d <- rep(NA_real_, n); d[bits==0] <- 2; d[bits==2] <- 1; d[bits==3] <- 0
  data.table(sample_id = sub("^N","", fam$V1), d = d)
}

snps <- c("rs9914422","rs744281","rs736866","rs7212001","rs726612","rs8070376","rs144461621")
dos <- list()
for (rs in snps) { s <- read_snp(rs); if (!is.null(s)) { setnames(s,"d",rs); dos[[rs]] <- s } }
present <- names(dos)
cat("SNPs present in CPGEA panel:", paste(present, collapse=", "), "\n")
cat("SNPs absent:", paste(setdiff(snps, present), collapse=", "), "\n")

M <- Reduce(function(a,b) merge(a,b,by="sample_id"), dos)
G <- as.matrix(M[, ..present])
R2 <- cor(G, use = "pairwise.complete.obs")^2
cat("\n=== pairwise r^2 ===\n"); print(round(R2,3))

df <- as.data.table(as.table(R2)); setnames(df, c("A","B","r2"))
df[, A := factor(A, levels = present)]
df[, B := factor(B, levels = rev(present))]
p <- ggplot(df, aes(A, B, fill = r2)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r2)), size = 3) +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0,1),
                      name = expression(r^2)) +
  labs(x = NULL, y = NULL,
       title = "KIF18B-locus haplotype LD (CPGEA)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(OUT, "FigS5_rs9914422_LD.pdf"), p, width = 6, height = 5.2, device = cairo_pdf)
ggsave(file.path(OUT, "FigS5_rs9914422_LD.png"), p, width = 6, height = 5.2, dpi = 300)
fwrite(as.data.table(R2, keep.rownames = "snp"), file.path(OUT, "FigS5_LD_matrix.csv"))
cat("\nSaved FigS5_rs9914422_LD.pdf\n")
