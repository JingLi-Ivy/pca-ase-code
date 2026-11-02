#!/usr/bin/env Rscript
# Does KIF18B carry prognostic information beyond proliferation?
# Cohorts: CPGEA, TCGA-PRAD (discovery) + MSKCC, DKFZ, GSE54460, GSE116918 (independent, as in 10b_external_validation.R).
# Proliferation: 7-gene score used throughout the manuscript (primary); Cuzick CCP 31-gene score (sensitivity; excludes KIF18B).
# Per cohort: Spearman(KIF18B, prolif); Cox BCR ~ KIF18B_z + prolif_z; delta C-index (clinical + prolif) -> (+ KIF18B), bootstrap.
suppressPackageStartupMessages({library(data.table); library(survival); library(jsonlite)})
E <- "/path/to/data_root/ASE/code/v8_round2/external"; setwd(E); set.seed(9)
P7 <- c("MKI67","PCNA","CCNB1","CCNB2","CDK1","CCNA2","AURKA")
CCP <- c("FOXM1","ASPM","TK1","PRC1","CDC20","BUB1B","PBK","DTL","CDKN3","RRM2","ASF1B","CEP55","CDK1","DLGAP5","SKA1","RAD51","KIF11",
         "BIRC5","RAD54L","CENPM","PCLAF","KIF20A","PTTG1","CDCA8","NUSAP1","PLK1","CDCA3","ORC6","CENPF","TOP2A","MCM10")
ALIAS <- c(PCLAF = "KIAA0101", ORC6 = "ORC6L", CDK1 = "CDC2")
G <- unique(c("KIF18B", P7, CCP))
sc <- function(M, genes) { g <- intersect(genes, colnames(M)); z <- scale(M[, g, drop = FALSE]); as.numeric(rowMeans(z, na.rm = TRUE)) }
fix <- function(M) { for (a in names(ALIAS)) if (!a %in% colnames(M) && ALIAS[[a]] %in% colnames(M)) colnames(M)[colnames(M) == ALIAS[[a]]] <- a; M }
gs <- function(x) { x <- as.character(x); s <- suppressWarnings(as.numeric(x)); i <- grepl("\\+", x); if (any(i)) s[i] <- vapply(strsplit(x[i], "+", fixed = TRUE), function(v) sum(as.numeric(v)), numeric(1)); s }
C <- list()

fp <- fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls"); tc <- grep("^T\\d+_WTS$", names(fp), value = TRUE)
f1 <- fp[gene_name %in% c(G, ALIAS)][!duplicated(gene_name)]
M <- t(log2(as.matrix(f1[, ..tc]) + 1)); colnames(M) <- f1$gene_name; M <- fix(M)
cl <- fread("/path/to/data_root/CPGEA/CPGEA_web/CPGEA_Clinical_Info.tsv")
d <- data.table(id = gsub("\\D", "", tc), K = M[, "KIF18B"], P7 = sc(M, P7), CCP = sc(M, CCP))
d <- merge(d, cl[, .(id = as.character(`Sequencing ID`), event = as.integer(toupper(BCR) == "YES"), time = suppressWarnings(as.numeric(sub("m$", "", `Month to surgery/BCR`))),
                     gleason = gs(GS), psa = suppressWarnings(as.numeric(PSA)))], by = "id")
C$CPGEA <- d[toupper(cl$BCR[match(id, cl$`Sequencing ID`)]) %in% c("YES","NO")]

e <- new.env(); load("/path/to/data_root/ASE/code/fig4_results/tcga_validation/TCGA_validation_workspace.RData", envir = e); cm <- e$clinical_matched
fu <- as.numeric(cm[["days_to_last_followup"]]); if (all(is.na(fu))) fu <- as.numeric(cm[["patient.days_to_last_followup"]])
bd <- as.numeric(cm[["patient.days_to_first_biochemical_recurrence"]]); st <- tolower(as.character(cm[["patient.biochemical_recurrence"]]))
ev <- ifelse(st %in% c("yes","true","1"), 1, ifelse(st %in% c("no","false","0"), 0, NA)); tv <- ifelse(ev == 1, dplyr::coalesce(bd, fu), fu)
psa_col <- grep("psa_value", names(cm), value = TRUE)[1]
gsn <- if ("GS_numeric" %in% names(cm)) cm$GS_numeric else cm[[grep("^gleason_score$", names(cm), value = TRUE)[1]]]
tb <- data.table(id = substr(cm$tcga_full_id, 1, 15), event = ev, time = tv / 30.44, gleason = suppressWarnings(as.numeric(gsn)),
                 psa = if (!is.na(psa_col)) suppressWarnings(as.numeric(cm[[psa_col]])) else NA_real_)
tp <- fread("/path/to/data_root/singlecell/scRNA_deconvolution/TCGA_PRAD_tumor.tpm", select = c(1, which(names(fread("/path/to/data_root/singlecell/scRNA_deconvolution/TCGA_PRAD_tumor.tpm", nrows = 0)) %in% c(G, ALIAS))))
setnames(tp, 1, "sid"); M <- log2(as.matrix(tp[, -1]) + 1); colnames(M) <- names(tp)[-1]; M <- fix(M)
tt <- data.table(id = substr(tp$sid, 1, 15), K = M[, "KIF18B"], P7 = sc(M, P7), CCP = sc(M, CCP))
tt <- tt[!duplicated(id)]
C$TCGA <- merge(tt, tb, by = "id")[!is.na(event)]
cat(sprintf("TCGA matched: %d samples (TPM ids like %s ; clinical ids like %s)\n", nrow(C$TCGA), tp$sid[1], cm$tcga_full_id[1]))

API <- "https://www.cbioportal.org/api"
post <- function(u, b) fromJSON(rawToChar(curl::curl_fetch_memory(paste0(API, u), handle = curl::handle_setheaders(curl::handle_setopt(curl::new_handle(),
  postfields = toJSON(b, auto_unbox = TRUE)), "Content-Type" = "application/json", Accept = "application/json"))$content))
ent <- post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", unique(c(G, ALIAS)))[, c("entrezGeneId","hugoGeneSymbol")]
for (x in list(c("MSKCC_Taylor2010","prad_mskcc","prad_mskcc_mrna_median_Zscores", "lin"), c("DKFZ_2018","prostate_dkfz_2018","prostate_dkfz_2018_rna_seq_mrna", "log"))) {
  md <- post(sprintf("/molecular-profiles/%s/molecular-data/fetch", x[3]), list(sampleListId = paste0(x[2], "_all"), entrezGeneIds = ent$entrezGeneId))
  md <- as.data.table(md)[, sym := ent$hugoGeneSymbol[match(entrezGeneId, ent$entrezGeneId)]]
  w <- dcast(md, sampleId + patientId ~ sym, value.var = "value", fun.aggregate = mean)
  M <- as.matrix(w[, -(1:2)]); if (x[4] == "log") M <- log2(M + 0.01); M <- fix(M)
  k <- data.table(sampleId = w$sampleId, patientId = w$patientId, K = M[, "KIF18B"], P7 = sc(M, P7), CCP = sc(M, CCP))
  c2 <- fread(sprintf("%s_clinical.csv", x[2]))
  if (x[1] == "MSKCC_Taylor2010") {
    z <- merge(k, c2, by = c("sampleId","patientId"))[SAMPLE_TYPE == "Primary" & DFS_STATUS != ""]
    C[[x[1]]] <- z[, .(id = patientId, K, P7, CCP, event = as.integer(DFS_STATUS == "1:Recurred"), time = as.numeric(DFS_MONTHS), gleason = as.numeric(GLEASON_SCORE), psa = NA_real_)]
  } else {
    z <- merge(k, c2, by = c("sampleId","patientId"))[INITIAL_TREATMENT == "RadP" & !is.na(BCR_STATUS)]
    C[[x[1]]] <- z[, .(K = mean(K), P7 = mean(P7), CCP = mean(CCP), event = as.integer(BCR_STATUS[1]), time = as.numeric(TIME_FROM_SURGERY_TO_BCR_LASTFU[1]),
                       gleason = gs(GLEASON_SCORE[1]), psa = as.numeric(PREOPERATIVE_PSA[1])), by = .(id = patientId)] } }

ser <- function(f) { L <- readLines(gzfile(f)); ch <- L[startsWith(L, "!Sample_characteristics_ch1")]
  M <- lapply(ch, function(l) gsub('"', "", strsplit(l, "\t")[[1]][-1])); key <- sapply(M, function(v) sub(":.*", "", v[1]))
  d <- as.data.table(setNames(lapply(M, function(v) sub("^[^:]*: ?", "", v)), make.names(key, unique = TRUE)))
  d[, gsm := gsub('"', "", strsplit(L[startsWith(L, "!Sample_geo_accession")], "\t")[[1]][-1])]; list(d = d, L = L) }
s <- ser("GSE54460_series_matrix.txt.gz")$d
f <- fread("GSE54460_FPKM-genes-TopHat2-106samples-12-4-13.txt.gz", header = FALSE)
ext <- unlist(f[2, -(1:2)]); f2 <- f[V2 %in% c(G, ALIAS)][!duplicated(V2)]
M <- t(log2(apply(as.matrix(f2[, -(1:2)]), 2, as.numeric) + 0.01)); colnames(M) <- f2$V2; M <- fix(M)
k <- data.table(ext = ext, K = M[, "KIF18B"], P7 = sc(M, P7), CCP = sc(M, CCP)); s[, ext := external.id]; s <- merge(s, k, by = "ext")
C$GSE54460_Long2014 <- s[, .(id = ext, K, P7, CCP, event = as.integer(bcr), time = fifelse(bcr == "1", suppressWarnings(as.numeric(months.to.bcr)), suppressWarnings(as.numeric(months.total.f.u))),
                             gleason = as.numeric(tgleason), psa = suppressWarnings(as.numeric(prepsa)))]

gp <- fread(cmd = "sed -n '/^ID\\t/,/^!platform_table_end/p' GPL25318_table.txt | grep -v '^!platform_table_end'")
gp <- gp[`Gene Symbol` %in% c(G, ALIAS) & grepl("^Sense \\(Fully Exonic\\)", Orientation)]
S2 <- ser("GSE116918_series_matrix.txt.gz"); L <- S2$L; b <- which(L == "!series_matrix_table_begin")
t2 <- fread(text = L[(b + 1):(which(L == "!series_matrix_table_end") - 1)]); setnames(t2, 1, "ID")
t2 <- merge(gp[, .(ID, sym = `Gene Symbol`)], t2, by = "ID")
agg <- t2[, lapply(.SD, mean), by = sym, .SDcols = setdiff(names(t2), c("ID","sym"))]
M <- t(as.matrix(agg[, -1])); colnames(M) <- agg$sym; M <- fix(M)
s2 <- S2$d; s2[, `:=`(K = M[gsm, "KIF18B"], P7 = sc(M, P7)[match(gsm, rownames(M))], CCP = sc(M, CCP)[match(gsm, rownames(M))])]
C$GSE116918_Jain2018 <- s2[, .(id = gsm, K, P7, CCP, event = as.integer(bcr.event..1.yes..0.no.), time = as.numeric(follow.up.time..bcr..months.),
                               gleason = as.numeric(gleason.grade), psa = suppressWarnings(as.numeric(psa..ng.ml.)))]

cidx <- function(f) summary(f)$concordance[1]
one <- function(nm, pr) { d <- copy(C[[nm]])[is.finite(K) & is.finite(get(pr)) & is.finite(time) & time > 0 & !is.na(event)]
  d[, `:=`(Kz = as.numeric(scale(K)), Pz = as.numeric(scale(get(pr))), g8 = as.integer(gleason >= 8))]
  usepsa <- sum(is.finite(d$psa) & d$psa > 0) >= 0.8 * nrow(d)
  u <- summary(coxph(Surv(time, event) ~ Kz, d))$coef; a <- summary(coxph(Surv(time, event) ~ Kz + Pz, d))$coef
  dd <- if (usepsa) d[is.finite(gleason) & is.finite(psa) & psa > 0] else d[is.finite(gleason)]
  f0 <- if (usepsa) Surv(time, event) ~ g8 + log(psa) + Pz else Surv(time, event) ~ g8 + Pz; f1 <- update(f0, . ~ . + Kz)
  ac <- summary(coxph(f1, dd))$coef["Kz", ]
  bs <- replicate(1000, { i <- sample(nrow(dd), replace = TRUE); x <- dd[i]; tryCatch(cidx(coxph(f1, x)) - cidx(coxph(f0, x)), error = function(e) NA) })
  data.table(cohort = nm, prolif = pr, n = nrow(d), events = sum(d$event), rho_K_prolif = cor(d$K, d[[pr]], method = "s"),
             HR_K = exp(u[1, 1]), P_K = u[1, 5], HR_K_adjProlif = exp(a["Kz", 1]), lo = exp(a["Kz", 1] - 1.96 * a["Kz", 3]), hi = exp(a["Kz", 1] + 1.96 * a["Kz", 3]),
             P_K_adjProlif = a["Kz", 5], HR_prolif = exp(a["Pz", 1]), P_prolif = a["Pz", 5], logHR = a["Kz", 1], se = a["Kz", 3],
             HR_K_adjClinProlif = exp(ac[1]), P_K_adjClinProlif = ac[5], dC = cidx(coxph(f1, dd)) - cidx(coxph(f0, dd)),
             dC_lo = quantile(bs, .025, na.rm = TRUE), dC_hi = quantile(bs, .975, na.rm = TRUE)) }
R <- rbindlist(lapply(c("P7","CCP"), function(pr) rbindlist(lapply(names(C), function(nm) tryCatch(one(nm, pr), error = function(e) { message(nm, ": ", conditionMessage(e)); NULL })))))
pool <- function(x) { w <- 1/x$se^2; fe <- sum(w*x$logHR)/sum(w); Q <- sum(w*(x$logHR-fe)^2); k <- nrow(x); t2 <- max(0, (Q-(k-1))/(sum(w)-sum(w^2)/sum(w)))
  wr <- 1/(x$se^2+t2); re <- sum(wr*x$logHR)/sum(wr); s <- sqrt(1/sum(wr)); sprintf("HR %.2f (%.2f-%.2f) P=%.3g I2=%.0f%%", exp(re), exp(re-1.96*s), exp(re+1.96*s), 2*pnorm(-abs(re/s)), 100*max(0,(Q-(k-1))/Q)) }
for (pr in c("P7","CCP")) { x <- R[prolif == pr]
  cat(sprintf("\n===== proliferation score: %s =====\n", pr))
  print(x[, .(cohort, n, events, rho = round(rho_K_prolif, 2), HR_K = round(HR_K, 2), P_K = signif(P_K, 2), HR_K_adjP = round(HR_K_adjProlif, 2),
              CI = sprintf("%.2f-%.2f", lo, hi), P_adjP = signif(P_K_adjProlif, 2), HR_prolif = round(HR_prolif, 2), P_prolif = signif(P_prolif, 2),
              HR_K_adjClinP = round(HR_K_adjClinProlif, 2), P_adjClinP = signif(P_K_adjClinProlif, 2), dC = sprintf("%+.3f [%+.3f,%+.3f]", dC, dC_lo, dC_hi))], width = 250)
  cat("Pooled KIF18B HR adjusted for proliferation, independent cohorts only:", pool(x[!cohort %in% c("CPGEA","TCGA")]), "\n")
  cat("Pooled KIF18B HR adjusted for proliferation, all six cohorts:", pool(x), "\n") }
fwrite(R, "kif18b_beyond_proliferation.csv")
