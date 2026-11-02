suppressPackageStartupMessages({library(data.table); library(ggplot2)})
OUT <- "/path/to/data_root/ASE/figures/pub_figs"
pca <- fread("/path/to/data_root/ASE/code/PCa_ASE_genes.csv")
fpkm <- fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls")
tc <- grep("^T\\d+_WTS$", names(fpkm), value=TRUE)
fpkm[, medFPKM := apply(.SD,1,median,na.rm=TRUE), .SDcols=tc]
m <- merge(pca[, .(gene_symbol, pct_sig_Tumor=as.numeric(pct_sig_Tumor), mean_ASE_Tumor=as.numeric(mean_ASE_Tumor))],
           fpkm[, .(gene_symbol=gene_name, medFPKM)], by="gene_symbol")
m <- m[is.finite(medFPKM) & medFPKM>0 & is.finite(pct_sig_Tumor)]
m[, expr_grp := cut(medFPKM, breaks=quantile(medFPKM, c(0,1/3,2/3,1)), include.lowest=TRUE,
                    labels=c("Low","Medium","High"))]
m <- m[!is.na(expr_grp)]
cat("n=",nrow(m),"\n"); print(m[, .(.N, med_pct=median(pct_sig_Tumor), med_ase=median(mean_ASE_Tumor)), by=expr_grp][order(expr_grp)])
kw <- kruskal.test(pct_sig_Tumor ~ expr_grp, data=m)
lab <- if(kw$p.value<2.2e-16) "P < 2.2e-16" else sprintf("P = %.2g", kw$p.value)
meds <- m[, .(med=median(pct_sig_Tumor)), by=expr_grp][order(expr_grp)]
medlab <- paste(sprintf("%.1f", meds$med), collapse=" / ")
ann <- sprintf("Kruskal-Wallis %s (medians %s%%)", lab, medlab)
cat(ann,"\n")
ytop <- quantile(m$pct_sig_Tumor, .99, na.rm=TRUE)
p <- ggplot(m, aes(expr_grp, pct_sig_Tumor, fill=expr_grp))+
  geom_boxplot(width=.6, outlier.size=.4)+
  scale_fill_manual(values=c(Low="#9ECAE1",Medium="#4292C6",High="#08519C"), guide="none")+
  annotate("text", x=2, y=ytop*1.08, label=ann, size=3.4)+
  coord_cartesian(ylim=c(0, ytop*1.12))+
  labs(x="Tumor expression tertile", y="Percent tumor samples with significant ASE")+
  theme_classic(base_size=12)
ggsave(file.path(OUT,"S4C_expression_stratification.pdf"), p, width=4.8, height=4.2)
ggsave(file.path(OUT,"S4C_expression_stratification.png"), p, width=4.8, height=4.2, dpi=300)
cat("DONE\n")
