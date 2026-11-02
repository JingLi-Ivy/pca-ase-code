#!/usr/bin/env Rscript
# =============================================================================
# Marker-based prostate lineage re-annotation + candidate-gene localization.
#
# The ism-pca deconvolution objects carry generic SingleR labels (Hepatocytes,
# Chondrocytes, iPS_cells...) that are implausible for prostate. Here we ignore
# those labels and assign each cell to a canonical prostate lineage by module
# scoring (Seurat::AddModuleScore, argmax over lineages), then localize the
# PCa-ASE candidate genes per lineage.
#
# Focus: is RBFOX3 epithelial (luminal/basal) or neural/stromal? Are the
# mitotic module (KIF18B/FOXM1/NEK2/TOP2A) and RNA-processing candidates
# (CELF2/EVI5/ITPR1/NCOA2) epithelial?
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(data.table); library(ggplot2) })

BASE   <- "/path/to/data_root/singlecell/scRNA_deconvolution"
OUTDIR <- "/path/to/data_root/ASE/figures/celltype_localization"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

GENES <- c("KIF18B","FOXM1","NEK2","TOP2A","RBFOX3","CELF2","EVI5","ITPR1","DRD1","NCOA2")

OBJS <- c(
  GSE176031 = "GSE176031/GSE176031_combined.rds",
  GSE181294 = "GSE181294/GSE181294_seurat_merged.rds",
  JGG       = "JGG/all.integrated_PCa.rds",
  GSE278936 = "GSE278936/combined_prostate_seurat.rds",
  GSE210358 = "GSE210358/GSE210358_seurat.rds"
)

MARK <- list(
  Luminal       = c("KLK3","KLK2","ACPP","NKX3-1","AR","KLK4","MSMB"),
  Basal         = c("KRT5","KRT14","TP63","KRT15","DST"),
  Club_Hillock  = c("SCGB1A1","SCGB3A1","LTF","PIGR","MMP7"),
  Endothelial   = c("VWF","PECAM1","CLDN5","FLT1","EGFL7"),
  Fibroblast    = c("DCN","LUM","COL1A1","COL1A2","PDGFRA"),
  SmoothMuscle  = c("ACTA2","MYH11","TAGLN","DES","CNN1"),
  Tcell         = c("CD3D","CD3E","CD8A","CD2","IL7R","CCL5"),
  Bcell_Plasma  = c("MS4A1","CD79A","MZB1","IGHG1","JCHAIN"),
  Myeloid       = c("LYZ","CD68","CD14","ITGAM","AIF1","C1QA"),
  Mast          = c("TPSAB1","CPA3","MS4A2"),
  Neuroendocrine= c("CHGA","CHGB","SYP","INSM1"),
  Neural        = c("PLP1","S100B","MPZ","NRXN1","SOX10")
)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), sprintf(...)))

all_long <- list(); comp_long <- list()
for (ds in names(OBJS)) {
  f <- file.path(BASE, OBJS[ds])
  if (!file.exists(f)) { log_msg("%s: missing %s; skip", ds, f); next }
  obj <- tryCatch(readRDS(f), error=function(e){log_msg("%s read fail: %s", ds, e$message); NULL})
  if (is.null(obj) || !inherits(obj,"Seurat")) { log_msg("%s: not seurat; skip", ds); next }
  res <- tryCatch({
    if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"
    obj <- tryCatch(JoinLayers(obj), error=function(e) obj)
    obj <- NormalizeData(obj, verbose=FALSE)
    sc_cols <- c()
    for (ln in names(MARK)) {
      g <- intersect(MARK[[ln]], rownames(obj))
      if (length(g) < 2) next
      ok2 <- tryCatch({
        obj <- AddModuleScore(obj, features=list(g), name=paste0("ms_",ln), seed=1, nbin=12, ctrl=50)
        TRUE
      }, error=function(e){ log_msg("%s: AddModuleScore %s failed (%s); skip lineage", ds, ln, conditionMessage(e)); FALSE })
      if (ok2) sc_cols[ln] <- paste0("ms_",ln,"1")
    }
    if (length(sc_cols) < 4) { log_msg("%s: too few lineage scores (%d); skip", ds, length(sc_cols)); return(NULL) }
    S <- as.matrix(obj@meta.data[, sc_cols, drop=FALSE]); colnames(S) <- names(sc_cols)
    best <- colnames(S)[max.col(S, ties.method="first")]
    best[apply(S,1,max) <= 0] <- "Unassigned"
    obj$lineage <- best
    tab <- as.data.table(table(lineage=best)); setnames(tab, "N", "n_cells"); tab[, dataset:=ds]
    comp_long[[ds]] <- tab
    log_msg("%s: re-annotated %d cells into %d lineages (%s)", ds, ncol(obj), length(unique(best)),
            paste(sort(unique(best)), collapse=","))
    g_present <- intersect(GENES, rownames(obj))
    expr <- tryCatch(GetAssayData(obj, layer="data"), error=function(e) GetAssayData(obj, slot="data"))
    expr <- as.matrix(expr[g_present, , drop=FALSE])
    for (g in g_present) {
      dt <- data.table(dataset=ds, gene=g, lineage=best, e=as.numeric(expr[g,]))
      all_long[[paste(ds,g,sep="|")]] <- dt[, .(mean_expr=mean(e), pct_expr=100*mean(e>0), n=.N),
                                            by=.(dataset,gene,lineage)]
    }
    TRUE
  }, error=function(e){ log_msg("%s: error (%s); skip", ds, conditionMessage(e)); NULL })
  rm(obj); gc()
}

if (length(all_long)==0) { log_msg("nothing produced; abort"); quit(status=1) }
L <- rbindlist(all_long); fwrite(L, file.path(OUTDIR,"reannotated_candidate_expression.csv"))
CC <- rbindlist(comp_long, fill=TRUE); fwrite(CC, file.path(OUTDIR,"reannotated_lineage_composition.csv"))
log_msg("wrote reannotated_candidate_expression.csv (%d rows)", nrow(L))

L[, epi := lineage %in% c("Luminal","Basal","Club_Hillock")]
epi_sum <- L[, .(epi_pct = weighted.mean(pct_expr[epi], n[epi]+1e-9),
                 nonepi_pct = weighted.mean(pct_expr[!epi], n[!epi]+1e-9)),
             by=.(gene)]
fwrite(epi_sum, file.path(OUTDIR,"reannotated_epi_vs_nonepi.csv"))

for (ds in unique(L$dataset)) {
  d <- L[dataset==ds]; d[, z:=scale(mean_expr)[,1], by=gene]
  p <- ggplot(d, aes(lineage, gene, size=pct_expr, color=z)) + geom_point() +
    scale_color_gradient2(low="#3B6EA8", mid="grey90", high="#B2182B", midpoint=0) +
    scale_size_continuous(range=c(0,7)) +
    labs(title=paste0("Marker-reannotated localization: ",ds), x=NULL, y=NULL,
         size="% expressing", color="scaled mean") +
    theme_bw(base_size=10) + theme(axis.text.x=element_text(angle=45,hjust=1))
  ggsave(file.path(OUTDIR, paste0("Fig_reannot_",ds,".pdf")), p,
         width=max(6,0.5*length(unique(d$lineage))+2), height=4.5)
}

log_msg("RBFOX3 epithelial vs non-epithelial (all datasets pooled):")
print(epi_sum[gene=="RBFOX3"])
log_msg("Done. Output: %s", OUTDIR)
