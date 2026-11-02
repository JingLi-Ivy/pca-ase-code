#!/usr/bin/env Rscript
# =============================================================================
# Cell-type localization of PCa-ASE candidate genes in annotated prostate
# cancer scRNA-seq datasets (reuses ism-pca annotated Seurat objects).
#
# Goal (revision P4): show whether candidate ASE genes are expressed in
# malignant epithelial cells vs stroma/immune/neuronal compartments, to rule
# out bulk-composition artifacts. Special focus on RBFOX3 (NeuN), which could
# be neuronal rather than epithelial.
#
# Output: per-dataset + combined per-cell-type mean expression and percent
# expressed for the candidate genes, plus a faceted dot-plot.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2)
})

BASE   <- "/path/to/data_root/singlecell/scRNA_deconvolution"
OUTDIR <- "/path/to/data_root/ASE/figures/celltype_localization"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

GENES <- c("KIF18B","FOXM1","NEK2","TOP2A","RBFOX3","CELF2","EVI5","ITPR1","DRD1","NCOA2")

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), sprintf(...)))

rds <- unique(c(
  Sys.glob(file.path(BASE, "*", "deconvolution_T", "*_seurat_annotated.rds")),
  Sys.glob(file.path(BASE, "*", "output", "*_annotated.rds")),
  file.path(BASE, "25965613/scRNA_total.rds")
))
rds <- rds[file.exists(rds)]
log_msg("Annotated Seurat rds present: %d", length(rds))
for (f in rds) cat("   ", f, "\n")
if (length(rds) == 0) { log_msg("No annotated rds found; abort."); quit(status=1) }

pick_celltype_col <- function(md) {
  cand <- c("cell_type","celltype","CellType","Cell_Type","cell.type",
            "annotation","Annotation","major_celltype","majortype",
            "cluster_annotation","seurat_annotations","cell_annotation")
  hit <- cand[cand %in% colnames(md)]
  if (length(hit)) return(hit[1])
  cc <- names(which(sapply(md, function(x) is.character(x) || is.factor(x))))
  cc <- setdiff(cc, c("orig.ident","sample","Sample","patient","Patient"))
  for (c in cc) { k <- length(unique(md[[c]])); if (k >= 3 && k <= 40) return(c) }
  return(NA_character_)
}

all_long <- list()
for (f in rds) {
  pdir <- basename(dirname(f))
  ds <- if (pdir %in% c("deconvolution_T","deconvolution_N","output")) basename(dirname(dirname(f))) else pdir
  if (pdir == "output") ds <- paste0(ds, "_int")
  obj <- tryCatch(readRDS(f), error=function(e){log_msg("read fail %s: %s", ds, e$message); NULL})
  if (is.null(obj)) next
  if (!inherits(obj, "Seurat")) { log_msg("%s: not a Seurat object (%s); skip", ds, class(obj)[1]); rm(obj); gc(); next }
  ok <- tryCatch({
    DefaultAssay(obj) <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
    md <- obj@meta.data
    ctcol <- pick_celltype_col(md)
    if (is.na(ctcol)) {
      if (length(levels(Idents(obj))) >= 3) { md$.__ct <- as.character(Idents(obj)); ctcol <- ".__ct" }
      else { log_msg("%s: no celltype column; skip", ds); ctcol <- NA_character_ }
    }
    if (!is.na(ctcol)) {
      log_msg("%s: celltype column = %s (%d types, %d cells)",
              ds, ctcol, length(unique(md[[ctcol]])), ncol(obj))
      g_present <- intersect(GENES, rownames(obj))
      if (length(g_present) == 0) { log_msg("%s: none of candidate genes present; skip", ds) }
      else {
        obj <- tryCatch(NormalizeData(obj, verbose=FALSE), error=function(e) obj)
        expr <- tryCatch(GetAssayData(obj, layer="data"), error=function(e) GetAssayData(obj, slot="data"))
        expr <- as.matrix(expr[g_present, , drop=FALSE])
        ct <- as.character(md[[ctcol]])
        for (g in g_present) {
          v <- as.numeric(expr[g, ])
          dt <- data.table(dataset=ds, gene=g, celltype=ct, e=v)
          agg <- dt[, .(mean_expr=mean(e), pct_expr=100*mean(e>0), n=.N), by=.(dataset,gene,celltype)]
          all_long[[paste(ds,g,sep="|")]] <- agg
        }
      }
    }
    TRUE
  }, error=function(e){ log_msg("%s: processing error (%s); skip", ds, conditionMessage(e)); FALSE })
  rm(obj); gc()
}

if (length(all_long) == 0) { log_msg("No expression extracted; abort."); quit(status=1) }
L <- rbindlist(all_long)
fwrite(L, file.path(OUTDIR, "celltype_candidate_expression.csv"))
log_msg("Wrote celltype_candidate_expression.csv (%d rows)", nrow(L))

for (ds in unique(L$dataset)) {
  d <- L[dataset==ds]
  d[, z := scale(mean_expr)[,1], by=gene]
  p <- ggplot(d, aes(x=celltype, y=gene, size=pct_expr, color=z)) +
    geom_point() +
    scale_color_gradient2(low="#3B6EA8", mid="grey90", high="#B2182B", midpoint=0) +
    scale_size_continuous(range=c(0,7)) +
    labs(title=paste0("Candidate ASE gene cell-type localization: ", ds),
         x=NULL, y=NULL, size="% expressing", color="scaled mean") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave(file.path(OUTDIR, paste0("Fig_celltype_", ds, ".pdf")), p,
         width=max(6, 0.4*length(unique(d$celltype))+3), height=4.5)
}

rb <- L[gene=="RBFOX3"]
if (nrow(rb)) {
  setorder(rb, dataset, -mean_expr)
  fwrite(rb, file.path(OUTDIR, "RBFOX3_celltype_summary.csv"))
  log_msg("RBFOX3 top cell types per dataset:")
  print(rb[, .SD[1:min(3,.N)], by=dataset][, .(dataset,celltype,mean_expr=round(mean_expr,3),pct_expr=round(pct_expr,1))])
}

log_msg("Done. Output: %s", OUTDIR)
