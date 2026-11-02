#!/usr/bin/env Rscript
# BayesPrism re-deconvolution of CPGEA tumour + normal bulk RNA-seq with the GSE185344 v2 single-cell reference,
# following the BayesPrism vignette (Chu et al. 2022 Nat Cancer) instead of a one-profile-per-type signature:
#   reference = raw single-cell counts (cells x genes); cell.type = celltype_major_v2 ("Epithelial" -> key "Epithelial_cells");
#   cell.state = celltype_refined_v2, with epithelial states further split by tissue of origin (tumour / benign sample) when the
#   metadata carry it; Unresolved removed; each state down-sampled to <= 2000 cells.
#   genes: cleanup.genes (Rb, Mrp, other_Rb, chrM, MALAT1, chrX, chrY), protein-coding only, then state markers
#   (get.exp.stat / select.marker, pval <= 0.01, lfc >= 0.1). mixture = raw bulk counts (read-only CPGEA matrices), samples x genes.
# (v2d: tumour/benign origin read from orig.ident, e.g. HYW_5742_Tumor / _Benign; v2c missed it because of a <= 6-level filter)
# Output: decon_v2c/CPGEA_{tumor,normal}_cell_fraction.tsv (same layout as decon_v2) + state-level fractions + QC comparisons.
suppressPackageStartupMessages({library(BayesPrism); library(Seurat); library(data.table)}); set.seed(32)
O <- Sys.getenv("DECON_OUT", "/path/to/data_root/ASE/code/v8_round2/decon_v2c"); dir.create(O, showWarnings = FALSE); setwd(O)
NC <- as.integer(Sys.getenv("NC", "16")); MAXC <- 2000L
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

seu <- readRDS("/path/to/data_root/singlecell/scRNA_deconvolution/TXNIP_v2/GSE185344/GSE185344_seurat_annotated_v2.rds")
m <- seu@meta.data; msg("cells %d ; metadata columns: %s", ncol(seu), paste(names(m), collapse = ", "))
for (cn in names(m)) if ((is.character(m[[cn]]) || is.factor(m[[cn]])) && length(unique(m[[cn]])) <= 30) { cat("==", cn, "\n"); print(sort(table(m[[cn]]), decreasing = TRUE)) }
type <- as.character(m$celltype_major_v2); type[type == "Epithelial"] <- "Epithelial_cells"
state <- if ("celltype_refined_v2" %in% names(m)) as.character(m$celltype_refined_v2) else type
oc <- NULL
for (cn in names(m)) { v <- tolower(as.character(m[[cn]])); u <- unique(v)
  if (length(u) >= 2 && length(u) <= 40 && any(grepl("tum|cancer|malig", u)) && any(grepl("benign|normal|adjacent", u))) { oc <- cn; break } }
if (!is.null(oc)) { org <- ifelse(grepl("benign|normal|adjacent", tolower(as.character(m[[oc]]))), "benign", "tumour")
  state[type == "Epithelial_cells"] <- paste(state[type == "Epithelial_cells"], org[type == "Epithelial_cells"], sep = "_")
  msg("tissue-of-origin column used for epithelial states: %s", oc) } else msg("no tumour/benign origin column found; epithelial states = refined labels only")
keep <- type != "Unresolved" & state != "Unresolved" & !is.na(type)
idx <- unlist(lapply(split(which(keep), state[keep]), function(i) if (length(i) > MAXC) sample(i, MAXC) else i))
msg("reference cells after down-sampling: %d", length(idx)); print(table(type[idx])); print(table(state[idx]))
ref <- t(as.matrix(GetAssayData(seu, assay = "RNA", layer = "counts")[, idx])); type <- type[idx]; state <- state[idx]; rm(seu); gc()

rd <- function(f) { d <- fread(f); g <- sub(".*\\|", "", d$gene_id); M <- as.matrix(d[, -1]); M <- rowsum(M, g); t(M) }
B <- rbind(rd("/path/to/data_root/raw/RNA_seq/CPGEA-RNA_seq/gene_count_matrix_tumor.csv"),
           rd("/path/to/data_root/raw/RNA_seq/CPGEA-RNA_seq/gene_count_matrix_normal.csv")[, , drop = FALSE])
msg("bulk samples %d (tumour %d, normal %d), genes %d", nrow(B), sum(grepl("^T", rownames(B))), sum(grepl("^N", rownames(B))), ncol(B))

ref <- cleanup.genes(input = ref, input.type = "count.matrix", species = "hs",
                     gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"), exp.cells = 5)
ref <- select.gene.type(ref, gene.type = "protein_coding")
st <- get.exp.stat(sc.dat = ref[, colSums(ref > 0) > 3], cell.type.labels = type, cell.state.labels = state,
                   pseudo.count = 0.1, cell.count.cutoff = 50, n.cores = NC)
ref <- select.marker(sc.dat = ref, stat = st, pval.max = 0.01, lfc.min = 0.1)
common <- intersect(colnames(ref), colnames(B)); msg("marker genes %d ; shared with bulk %d", ncol(ref), length(common))

pr <- new.prism(reference = ref[, common], mixture = B[, common], input.type = "count.matrix",
                cell.type.labels = type, cell.state.labels = state, key = "Epithelial_cells", outlier.cut = 0.01, outlier.fraction = 0.1)
res <- run.prism(prism = pr, n.cores = NC)
th <- get.fraction(bp = res, which.theta = "final", state.or.type = "type")
ths <- get.fraction(bp = res, which.theta = "first", state.or.type = "state")
saveRDS(res, file.path(O, "bayesprism_v2c_result.rds"))
wr <- function(M, f) fwrite(data.table(sample = rownames(M), M), f, sep = "\t")
wr(th[grepl("^T", rownames(th)), , drop = FALSE], "CPGEA_tumor_cell_fraction.tsv"); wr(th[grepl("^N", rownames(th)), , drop = FALSE], "CPGEA_normal_cell_fraction.tsv")
wr(ths, "CPGEA_state_fraction.tsv")

for (ti in c("T", "N")) { x <- th[grepl(paste0("^", ti), rownames(th)), , drop = FALSE]; msg("[%s] median fractions:", ti); print(round(apply(x, 2, median), 3)) }
old <- rbind(fread("/path/to/data_root/ASE/code/v8_round2/decon_v2/CPGEA_tumor_cell_fraction.tsv"), fread("/path/to/data_root/ASE/code/v8_round2/decon_v2/CPGEA_normal_cell_fraction.tsv"))
setnames(old, 1, "sample"); E <- log2(B[, intersect(c("EPCAM", "KRT8", "KRT18", "CDH1", "CLDN4", "KRT19"), colnames(B))] / rowSums(B) * 1e6 + 1)
mk <- rowMeans(scale(E))
for (ti in c("T", "N")) { s <- rownames(th)[grepl(paste0("^", ti), rownames(th))]
  msg("[%s] epithelial: Spearman with previous signature run %.2f ; with marker score: new %.2f vs previous %.2f", ti,
      cor(th[s, "Epithelial_cells"], old[match(s, sample), Epithelial_cells], method = "s", use = "c"),
      cor(th[s, "Epithelial_cells"], mk[s], method = "s"), cor(old[match(s, sample), Epithelial_cells], mk[s], method = "s", use = "c")) }
if (any(grepl("tumour", colnames(ths)))) { tc <- grep("_tumour$", colnames(ths), value = TRUE)
  msg("tumour-derived epithelial state fraction, median tumour %.3f vs normal %.3f", median(rowSums(ths[grepl("^T", rownames(ths)), tc, drop = FALSE])), median(rowSums(ths[grepl("^N", rownames(ths)), tc, drop = FALSE]))) }
msg("DONE")
