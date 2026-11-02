#!/usr/bin/env Rscript
# 37b: rerun of 37 with (i) removal of overflowed counts (negative or depth >= 32767; 16-bit overflow in upstream count files),
# (ii) one output file per patient and mode (resumable, progress visible), (iii) pass 1 = one-sample mode for all patients,
# pass 2 = two-sample mode (slow). Original header of 37:
# Independent re-test of Delta-ASE with MBASED (Mayba et al. 2014 Genome Biol; Bioconductor), gene-level ASE with internal
# pseudo-phasing of the gene's SNPs (isPhased = FALSE), per patient.
#   (1) one-sample mode, run separately on tumour and normal counts of the same SNPs -> gene major-allele frequency (MAF) in each
#       tissue; symmetric comparison of MAF_T vs MAF_N (no orientation chosen by either tissue);
#   (2) two-sample mode (MBASED's tumour-vs-normal differential ASE), run in both orders (tumour first / normal first): the
#       sample-1 major allele orients the difference, so each order favours its own sample-1; comparing the two runs cancels this.
# Sites: paired balanced-diploid sites (18) after the mappability filters of 27; per patient x gene, SNPs >= 300 bp apart, <= 20 deepest; numSim = 0.
suppressPackageStartupMessages({library(data.table); library(MBASED); library(SummarizedExperiment); library(GenomicRanges); library(BiocParallel); library(parallel)})
V <- "/path/to/data_root/ASE/code/v8_round2"; D <- file.path(V, "delta_ase"); F <- file.path(V, "mapfilter"); O <- file.path(V, "delta_mbased"); dir.create(O, showWarnings = FALSE)
NC <- as.integer(Sys.getenv("NC", "16")); msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

B <- fread(file.path(D, "paired_balanced_sites.csv.gz"))
B[, `:=`(snp = paste(chrom, pos, ref, alt, sep = ":"), chrom = as.character(chrom), eg = sub("\\..*", "", eg))]
sd <- fread(file.path(F, "genomicSuperDups.txt.gz"), select = 2:4, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
bl <- fread(file.path(F, "hg38-blacklist.v2.bed.gz"), select = 1:3, col.names = c("chr", "s", "e"))[, chrom := sub("^chr", "", chr)]
ov <- function(reg) { q <- unique(B[, .(chrom, s = pos, e = pos)]); r <- unique(reg[, .(chrom, s = s + 1L, e)]); setkey(r, chrom, s, e)
  foverlaps(q, r, nomatch = 0L)[, unique(paste(chrom, i.s))] }
bad <- c(ov(sd), ov(bl))
gtf <- fread(cmd = "zcat /path/to/data_root/reference/hg38/Homo_sapiens.GRCh38.110.gtf.gz | awk -F'\\t' '$3==\"gene\"' | cut -f9", header = FALSE, sep = "\t")
pseudo <- gtf[grepl('gene_biotype "[a-z_]*pseudogene"', V1), sub('.*gene_id "([^"]+)".*', "\\1", V1)]
mono <- unique(B, by = c("pid", "snp"))[, .(devN = abs(sum(rnN) / sum(rnN + anN) - 0.5)), by = snp][devN > 0.35, snp]
B <- B[!(paste(chrom, pos) %chin% bad | (chrom == "6" & pos >= 28510120 & pos <= 33480577) | eg %chin% pseudo | snp %chin% mono)]
B <- unique(B, by = c("pid", "snp")); B[, pid := as.character(pid)]
nb <- nrow(B); B <- B[rnT >= 0 & anT >= 0 & rnN >= 0 & anN >= 0 & rnT + anT < 32767 & rnN + anN < 32767]
msg("removed %d rows with overflowed counts", nb - nrow(B))
MAXSNP <- 20L; B[, dep := rnT + anT + rnN + anN]; setorder(B, pid, eg, -dep)
B <- B[, { kp <- logical(.N); kept <- integer(0); for (j in seq_len(.N)) if (length(kept) < MAXSNP && (!length(kept) || all(abs(pos[j] - kept) >= 300L))) { kp[j] <- TRUE; kept <- c(kept, pos[j]) }; .SD[kp] }, by = .(pid, eg)]
B[, `:=`(rnT = as.integer(rnT), anT = as.integer(anT), rnN = as.integer(rnN), anN = as.integer(anN))]
msg("patient-SNP observations %s ; patients %d", format(nrow(B), big.mark = ","), uniqueN(B$pid))

mk_se <- function(d, a1, a2) {
  gr <- GRanges(paste0("chr", d$chrom), IRanges(d$pos, width = 1), aseID = d$eg, allele1 = d$ref, allele2 = d$alt); names(gr) <- paste(d$snp, d$eg)
  a1 <- as.matrix(a1); a2 <- as.matrix(a2); storage.mode(a1) <- "integer"; storage.mode(a2) <- "integer"
  SummarizedExperiment(assays = list(lociAllele1Counts = a1, lociAllele2Counts = a2), rowRanges = gr) }
run1 <- function(d, a1, a2) { x <- runMBASED(ASESummarizedExperiment = mk_se(d, a1, a2), isPhased = FALSE, numSim = 0, BPPARAM = SerialParam())
  data.table(eg = rownames(assays(x)$majorAlleleFrequency), maf = assays(x)$majorAlleleFrequency[, 1], p = assays(x)$pValueASE[, 1]) }
run2 <- function(d, a1, a2) { x <- runMBASED(ASESummarizedExperiment = mk_se(d, a1, a2), isPhased = FALSE, numSim = 0, BPPARAM = SerialParam())
  data.table(eg = rownames(assays(x)$majorAlleleFrequencyDifference), diff = assays(x)$majorAlleleFrequencyDifference[, 1], p = assays(x)$pValueASE[, 1]) }
dir.create(file.path(O, "one"), showWarnings = FALSE); dir.create(file.path(O, "two"), showWarnings = FALSE)
pass <- function(mode) { pids <- unique(B$pid)
  invisible(mclapply(pids, function(p) { f <- file.path(O, mode, paste0(p, ".csv")); if (file.exists(f)) return(NULL)
    d <- B[pid == p]; if (nrow(d) < 5) return(NULL); t0 <- Sys.time()
    r <- tryCatch(if (mode == "one") merge(run1(d, matrix(d$rnT), matrix(d$anT)), run1(d, matrix(d$rnN), matrix(d$anN)), by = "eg", suffixes = c("_T", "_N"))
                  else merge(run2(d, cbind(d$rnT, d$rnN), cbind(d$anT, d$anN)), run2(d, cbind(d$rnN, d$rnT), cbind(d$anN, d$anT)), by = "eg", suffixes = c("_Tfirst", "_Nfirst")),
                  error = function(e) { message(p, " ", mode, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) fwrite(r[, pid := p], f)
    message(sprintf("%s pass %s patient %s: %d genes, %.0f s", format(Sys.time(), "%H:%M"), mode, p, if (is.null(r)) 0L else nrow(r), as.numeric(Sys.time() - t0, units = "secs"))) }, mc.cores = NC)) }
summ <- function() {
  X <- rbindlist(lapply(list.files(file.path(O, "one"), full.names = TRUE), fread))[is.finite(maf_T) & is.finite(maf_N)]
  pp <- X[, .(d = mean(maf_T - maf_N)), by = pid]
  msg("\n(1) one-sample, %d patients, %s patient-genes: MAF tumour %.4f vs normal %.4f ; genes MAF_T > MAF_N %.1f%% ; patients mean MAF_T > MAF_N %.0f%% (Wilcoxon P=%.2g)",
      uniqueN(X$pid), format(nrow(X), big.mark = ","), X[, mean(maf_T)], X[, mean(maf_N)], 100 * X[, mean(maf_T > maf_N)], 100 * pp[, mean(d > 0)], wilcox.test(pp$d)$p.value)
  for (th in c(0.05, 0.01)) { a <- X[p_T < th | p_N < th]; msg("    genes with ASE (P<%g) in either tissue: n=%s, MAF_T > MAF_N %.1f%%", th, format(nrow(a), big.mark = ","), 100 * a[, mean(maf_T > maf_N)]) }
  fl <- list.files(file.path(O, "two"), full.names = TRUE)
  if (length(fl)) { R <- rbindlist(lapply(fl, fread)); msg("(2) two-sample, %d patients:", uniqueN(R$pid))
    for (th in c(0.05, 0.01, 0.001)) { st <- R[p_Tfirst < th, .N]; sn <- R[p_Nfirst < th, .N]
      msg("    differential ASE P<%g: tumour-first %s genes ; normal-first %s genes ; ratio %.2f", th, format(st, big.mark = ","), format(sn, big.mark = ","), st / sn) } } }
pass("one"); summ(); msg("pass one DONE")
pass("two"); summ(); msg("DONE")
