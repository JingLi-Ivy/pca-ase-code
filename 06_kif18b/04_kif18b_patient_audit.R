#!/usr/bin/env Rscript
# KIF18B patient evidence table: one row per patient with tumor RNA.
suppressPackageStartupMessages(library(data.table))
O <- "/path/to/data_root/ASE/code/v8_round2/kif18b_audit"; setwd(O)
GID <- "ENSG00000186185"
RS_POS <- 44947872
pid_of <- function(x) gsub("\\D", "", x)

g <- fread("rs9914422.raw"); gc <- grep("rs9914422", names(g), value = TRUE)
geno <- g[, .(pid = pid_of(FID), tissue = substr(FID, 1, 1), cnt = get(gc), counted = sub("rs9914422_", "", gc))]
if (all(geno$tissue %in% c("0","1","2","3","4","5","6","7","8","9"))) geno[, tissue := "?"]
geno <- geno[order(pid, tissue)][, .SD[1], by = pid]
cat(sprintf("genotyped patients: %d (allele counted: %s)\n", nrow(geno), geno$counted[1]))

fp <- fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls")
fp[, gene_id := sub("\\..*", "", gene_id)]
tc <- grep("^T\\d+_WTS$", names(fp), value = TRUE)
getg <- function(col, val) as.numeric(fp[get(col) == val][1, ..tc])
prolif <- c("MKI67","PCNA","CCNB1","CCNB2","CDK1","CCNA2","AURKA")
pm <- log2(as.matrix(fp[gene_name %in% prolif, ..tc]) + 1)
ex <- data.table(pid = pid_of(tc), FOXM1 = getg("gene_name","FOXM1"), KIF18B = getg("gene_id",GID),
                 prolif_score = as.numeric(scale(colMeans(t(scale(t(pm)))))))
ex[, FOXM1_tertile := as.integer(cut(FOXM1, c(-Inf, quantile(FOXM1, c(1/3,2/3)), Inf)))]

fs <- list.files("/path/to/data_root/ASE", pattern = "^T\\d+_ase_results\\.csv$", full.names = TRUE)
sen <- rbindlist(lapply(fs, function(f) { d <- fread(f)[gene_id == GID, .(a, pval, fdr)]; if (!nrow(d)) return(NULL); d[, pid := pid_of(basename(f))]; d }))
ag <- fread("/path/to/data_root/ASE/code/ase_gene_sample_agg.csv.gz")[gene_id == GID, .(pid = pid_of(sample_id), allelic_depth = dep)]

r <- fread("kif18b_rna_snps.csv")[gene_id == GID]
r[, pid := pid_of(sample_id)]
rs <- r[, .(rna_n_snp = .N, rna_total_reads = sum(N),
            rna_snps = paste(sprintf("%d:%s>%s(%d/%d)", start, ref, alt, ref.matches, N), collapse = ";"),
            rna_rs9914422_het = any(start == RS_POS), n_sig_snp = sum(ase_status == "Significant_ASE")), by = pid]

p <- fread("kif18b_dna_pileup.csv", header = FALSE,
           col.names = c("sample","chr","pos","ref","alt","NR","NA_","NE","ND","TR","TA","TE","TD"))
p[, pid := pid_of(sample)]
p <- p[ref != "." & (NR + NA_) >= 15][, nvaf := NA_ / (NR + NA_)][nvaf >= 0.25 & nvaf <= 0.75]
p[, tvaf := TA / pmax(TR + TA, 1)]
dna <- p[, .(dna_n_het = .N, dna_tumor_depth = sum(TR + TA),
             dna_tumor_imb = round(median(abs(tvaf - 0.5)), 3), dna_normal_imb = round(median(abs(nvaf - 0.5)), 3),
             dna_rs9914422 = if (any(pos == RS_POS)) sprintf("N %d/%d T %d/%d", NA_[pos==RS_POS], NR[pos==RS_POS]+NA_[pos==RS_POS], TA[pos==RS_POS], TR[pos==RS_POS]+TA[pos==RS_POS]) else NA_character_),
         by = pid]
m <- merge(r[, .(pid, pos = start, rna_ref = ref.matches, rna_n = N)], p[, .(pid, pos, TR, TA)], by = c("pid","pos"))
dr <- m[rna_n >= 1, .(n_matched_sites = .N,
                      rna_ref_frac = round(sum(rna_ref) / sum(rna_n), 3),
                      dna_ref_frac = round(sum(TR) / sum(TR + TA), 3)), by = pid]

fac <- fread("/path/to/data_root/variants/batch_43/facets/annotated_cnv_output.tsv")[gene_name == "KIF18B"]
fac <- fac[order(sample_id, -num_mark)][, .SD[1], by = sample_id][, .(pid = pid_of(sample_id), tcn, lcn, cf, cnv_type)]
pur <- fread("/path/to/data_root/variants/batch_43/facets/Table_S1_PRAD.txt")[, .(pid = as.character(SAMPLE_ID), purity = round(as.numeric(purity), 3))]
dec <- fread("/path/to/data_root/singlecell/scRNA_deconvolution/GSE278936/deconvolution_T/output/CPGEA_cell_fraction.tsv")
idc <- names(dec)[1]; dec <- dec[, .(pid = pid_of(get(idc)), epithelial = round(Epithelial_cells, 3))]

cl <- fread("/path/to/data_root/CPGEA/CPGEA_web/CPGEA_Clinical_Info.tsv")
cl <- cl[, .(pid = as.character(`Sequencing ID`), GS, pT = `Pathology T Stage`, BCR, BCR_months = `Month to surgery/BCR`)]

dtf <- "/path/to/data_root/ASE/code/v8_round2/kif18b_dt/dt_counts.tsv"
dt <- if (file.exists(dtf)) fread(dtf, colClasses = c(pid = "character"))[, .(pid, KIF18B_DT = DT_reads / mapped * 1e6, KIF18B_rpm = KIF18B_reads / mapped * 1e6)] else data.table(pid = character(), KIF18B_DT = numeric(), KIF18B_rpm = numeric())
tab <- Reduce(function(a, b) merge(a, b, by = "pid", all.x = TRUE),
              list(ex, dt, geno[, .(pid, rs9914422_dosage = cnt)], sen, ag, rs, dna, dr, fac, pur, dec, cl))
tab[, rs9914422_class := fcase(rs9914422_dosage == 1, "het", rs9914422_dosage %in% c(0, 2), "hom", default = NA_character_)]
tab[, sen_evaluable := !is.na(fdr)][, sig_ASE := !is.na(fdr) & fdr < 0.10]
setcolorder(tab, c("pid","sen_evaluable","sig_ASE","a","fdr","allelic_depth","rs9914422_dosage","rs9914422_class",
                   "FOXM1","FOXM1_tertile","prolif_score","KIF18B","KIF18B_DT","tcn","lcn","cnv_type","cf","purity","epithelial"))
fwrite(tab[order(-sig_ASE, -sen_evaluable, FOXM1_tertile)], "kif18b_patient_audit_all.csv")
cat(sprintf("\nrows=%d  Sen evaluable=%d  sig=%d  genotyped among evaluable=%d\n", nrow(tab), sum(tab$sen_evaluable), sum(tab$sig_ASE), tab[sen_evaluable & !is.na(rs9914422_class), .N]))

cat("\n===== 9 KIF18B ASE carriers + T67 =====\n")
print(tab[sig_ASE | pid == "67", .(pid, fdr = signif(fdr, 2), a = round(a, 2), dep = allelic_depth, geno = rs9914422_class, FOXM1 = round(FOXM1, 1), tert = FOXM1_tertile,
         tcn, lcn, cnv_type, purity, epi = epithelial, KIF18B = round(KIF18B, 2), DT = round(KIF18B_DT, 2), rna_ref_frac, dna_ref_frac, n_matched_sites, GS, BCR)])

cat("\n===== Local copy state: carriers vs non-carriers (evaluable) =====\n")
tab[, cn_state := fcase(is.na(tcn), "no_call", tcn == 2 & lcn == 1, "balanced", tcn == 2 & lcn == 0, "CN-LOH", tcn < 2, "loss", tcn > 2, "gain", default = "other")]
print(dcast(tab[sen_evaluable == TRUE, .N, by = .(sig_ASE, cn_state)], cn_state ~ sig_ASE, value.var = "N", fill = 0))
cat("Fisher (imbalanced vs balanced, carriers vs non): p=",
    signif(fisher.test(table(tab[sen_evaluable == TRUE & cn_state != "no_call", .(sig_ASE, cn_state != "balanced")]))$p.value, 3), "\n")

cat("\n===== KIF18B vs KIF18B-DT =====\n")
cat(sprintf("KIF18B-DT median FPKM=%.2f; Spearman(KIF18B, DT)=%.3f p=%.2g; (FOXM1, DT)=%.3f p=%.2g\n",
    median(tab$KIF18B_DT, na.rm = TRUE),
    cor(tab$KIF18B, tab$KIF18B_DT, method = "s", use = "c"), cor.test(tab$KIF18B, tab$KIF18B_DT, method = "s", exact = FALSE)$p.value,
    cor(tab$FOXM1, tab$KIF18B_DT, method = "s", use = "c"), cor.test(tab$FOXM1, tab$KIF18B_DT, method = "s", exact = FALSE)$p.value))
for (y in c("KIF18B","KIF18B_DT")) {
  z <- tab[!is.na(rs9914422_dosage)]
  f <- lm(log2(get(y) + 0.01) ~ rs9914422_dosage, data = z); s <- summary(f)$coef[2, ]
  cat(sprintf("%s ~ rs9914422 (counted %s), n=%d: beta=%.3f p=%.3g ; medians by dosage 0/1/2 = %s\n", y, geno$counted[1], nrow(z), s[1], s[4],
      paste(round(z[, median(get(y)), by = rs9914422_dosage][order(rs9914422_dosage)]$V1, 2), collapse = "/")))
}
cat("\n===== T67 relative position =====\n")
for (v in c("FOXM1","prolif_score","KIF18B","KIF18B_DT","purity","epithelial")) {
  x <- tab[[v]]; cat(sprintf("%-12s T67=%.3f  rank %d/%d  z=%.2f\n", v, tab[pid=="67"][[v]], rank(-x, na.last="keep")[tab$pid=="67"], sum(!is.na(x)), (tab[pid=="67"][[v]] - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)))
}
