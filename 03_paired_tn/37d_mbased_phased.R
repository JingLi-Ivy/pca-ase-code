#!/usr/bin/env Rscript
# 37d: MBASED with real phase (isPhased = TRUE) instead of pseudo-phasing. Sites, filters and per patient x gene SNP selection as 37b,
# restricted to patient-SNPs that are phased heterozygous (0|1 or 1|0) in the Beagle 5.5 + 1000G genome-wide phasing (31_phase_genome.sh);
# SNP selection (>= 300 bp apart, <= MAXSNP deepest) is done after this restriction. Allele 1 = haplotype 1 allele.
# Purpose: show that the 37b result does not depend on pseudo-phasing (which inflates MAF and lets sample 1 choose the major allele).
#   (1) one-sample, tumour and normal separately: MAF_T vs MAF_N per patient (as 37b); haplotype-1 frequency (signed) in each tissue;
#       among genes with ASE (P<0.05) in both tissues, how often the same haplotype is major in tumour and normal.
#   (2) two-sample, all genes, both orders (MBASED still takes the major allele from sample 1, phased or not).
# MBASED does not model overdispersion: counts of significant genes are inflated; read directions and patient-level comparisons only.
#   (2) two-sample mode (MBASED's tumour-vs-normal differential ASE), run in both orders (tumour first / normal first): the
#       sample-1 major allele orients the difference, so each order favours its own sample-1; comparing the two runs cancels this.
# Sites: paired balanced-diploid sites (18) after the mappability filters of 27; per patient x gene, SNPs >= 300 bp apart, <= 20 deepest; numSim = 0.
suppressPackageStartupMessages({library(data.table); library(MBASED); library(SummarizedExperiment); library(GenomicRanges); library(BiocParallel); library(parallel)})
V <- "/path/to/data_root/ASE/code/v8_round2"; D <- file.path(V, "delta_ase"); F <- file.path(V, "mapfilter"); O <- file.path(V, "delta_mbased_phased"); dir.create(O, showWarnings = FALSE)
PH <- file.path(V, "phase_genome"); MAXSNP <- as.integer(Sys.getenv("MAXSNP", "20"))
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
pf <- file.path(O, "positions.tsv"); fwrite(unique(B[, .(chrom, pos)])[order(as.integer(chrom), pos)], pf, sep = "\t", col.names = FALSE)
GTL <- rbindlist(lapply(1:22, function(c) { f <- file.path(PH, sprintf("phased_chr%d.vcf.gz", c)); if (!file.exists(f)) return(NULL)
  smp <- system(sprintf("bcftools query -l %s", f), intern = TRUE)
  x <- fread(cmd = sprintf("bcftools query -R %s -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' %s", pf, f), header = FALSE, sep = "\t", colClasses = "character")
  if (!nrow(x)) return(NULL); setnames(x, c("chrom", "pos", "vref", "valt", gsub("\\D", "", smp)))
  melt(x, id.vars = c("chrom", "pos", "vref", "valt"), variable.name = "pid", value.name = "gt")[gt %in% c("0|1", "1|0")] }))
GTL[, `:=`(pos = as.integer(pos), pid = as.character(pid), h1 = substr(gt, 1, 1))]
n0 <- nrow(B); B <- merge(B, GTL[, .(chrom, pos, pid, vref, valt, h1)], by = c("chrom", "pos", "pid"))
B[, h1allele := fifelse(h1 == "0", vref, valt)]
B <- B[(vref == ref & valt == alt) | (vref == alt & valt == ref)]
B[, `:=`(a1T = fifelse(h1allele == ref, rnT, anT), a2T = fifelse(h1allele == ref, anT, rnT),
         a1N = fifelse(h1allele == ref, rnN, anN), a2N = fifelse(h1allele == ref, anN, rnN))]
msg("patient-SNPs phased heterozygous with matching alleles: %s of %s (%.1f%%)", format(nrow(B), big.mark = ","), format(n0, big.mark = ","), 100 * nrow(B) / n0)
B[, dep := rnT + anT + rnN + anN]; setorder(B, pid, eg, -dep)
B <- B[, { kp <- logical(.N); kept <- integer(0); for (j in seq_len(.N)) if (length(kept) < MAXSNP && (!length(kept) || all(abs(pos[j] - kept) >= 300L))) { kp[j] <- TRUE; kept <- c(kept, pos[j]) }; .SD[kp] }, by = .(pid, eg)]
B[, `:=`(rnT = as.integer(rnT), anT = as.integer(anT), rnN = as.integer(rnN), anN = as.integer(anN))]
msg("patient-SNP observations %s ; patients %d", format(nrow(B), big.mark = ","), uniqueN(B$pid))
B[, `:=`(a1T = as.integer(a1T), a2T = as.integer(a2T), a1N = as.integer(a1N), a2N = as.integer(a2N))]

mk_se <- function(d, a1, a2) {
  gr <- GRanges(paste0("chr", d$chrom), IRanges(d$pos, width = 1), aseID = d$eg, allele1 = d$h1allele, allele2 = fifelse(d$h1allele == d$ref, d$alt, d$ref)); names(gr) <- paste(d$snp, d$eg)
  a1 <- as.matrix(a1); a2 <- as.matrix(a2); storage.mode(a1) <- "integer"; storage.mode(a2) <- "integer"
  SummarizedExperiment(assays = list(lociAllele1Counts = a1, lociAllele2Counts = a2), rowRanges = gr) }
a1maj <- function(x) { l <- metadata(x)$locusSpecificResults; v <- assays(l)$allele1IsMajor[, 1]
  tapply(v, as.character(mcols(rowRanges(l))$aseID), function(z) z[1]) }
run1 <- function(d, a1, a2) { x <- runMBASED(ASESummarizedExperiment = mk_se(d, a1, a2), isPhased = TRUE, numSim = 0, BPPARAM = SerialParam())
  r <- data.table(eg = rownames(assays(x)$majorAlleleFrequency), maf = assays(x)$majorAlleleFrequency[, 1], p = assays(x)$pValueASE[, 1])
  m <- a1maj(x); r[, h1f := fifelse(m[eg], maf, 1 - maf)]; r }
run2 <- function(d, a1, a2) { x <- runMBASED(ASESummarizedExperiment = mk_se(d, a1, a2), isPhased = TRUE, numSim = 0, BPPARAM = SerialParam())
  data.table(eg = rownames(assays(x)$majorAlleleFrequencyDifference), diff = assays(x)$majorAlleleFrequencyDifference[, 1], p = assays(x)$pValueASE[, 1]) }
dir.create(file.path(O, "one"), showWarnings = FALSE); dir.create(file.path(O, "two"), showWarnings = FALSE)
pass <- function(mode) { pids <- unique(B$pid)
  invisible(mclapply(pids, function(p) { f <- file.path(O, mode, paste0(p, ".csv")); if (file.exists(f)) return(NULL)
    d <- B[pid == p]; if (nrow(d) < 5) return(NULL); t0 <- Sys.time()
    r <- tryCatch(if (mode == "one") merge(run1(d, matrix(d$a1T), matrix(d$a2T)), run1(d, matrix(d$a1N), matrix(d$a2N)), by = "eg", suffixes = c("_T", "_N"))
                  else merge(run2(d, cbind(d$a1T, d$a1N), cbind(d$a2T, d$a2N)), run2(d, cbind(d$a1N, d$a1T), cbind(d$a2N, d$a2T)), by = "eg", suffixes = c("_Tfirst", "_Nfirst")),
                  error = function(e) { message(p, " ", mode, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) fwrite(r[, pid := p], f)
    message(sprintf("%s pass %s patient %s: %d genes, %.0f s", format(Sys.time(), "%H:%M"), mode, p, if (is.null(r)) 0L else nrow(r), as.numeric(Sys.time() - t0, units = "secs"))) }, mc.cores = NC)) }
summ <- function() {
  X <- rbindlist(lapply(list.files(file.path(O, "one"), full.names = TRUE), fread))[is.finite(maf_T) & is.finite(maf_N)]
  pp <- X[, .(d = mean(maf_T - maf_N)), by = pid]
  msg("\n(1) one-sample (phased), %d patients, %s patient-genes: MAF tumour %.4f vs normal %.4f ; genes MAF_T > MAF_N %.1f%% ; patients mean MAF_T > MAF_N %.0f%% (Wilcoxon P=%.2g)",
      uniqueN(X$pid), format(nrow(X), big.mark = ","), X[, mean(maf_T)], X[, mean(maf_N)], 100 * X[, mean(maf_T > maf_N)], 100 * pp[, mean(d > 0)], wilcox.test(pp$d)$p.value)
  for (th in c(0.05, 0.01)) { a <- X[p_T < th | p_N < th]; msg("    genes with ASE (P<%g) in either tissue: n=%s, MAF_T > MAF_N %.1f%%", th, format(nrow(a), big.mark = ","), 100 * a[, mean(maf_T > maf_N)]) }
  b <- X[p_T < 0.05 & p_N < 0.05]; msg("    genes with ASE (P<0.05) in both tissues: n=%s, same major haplotype in T and N %.1f%%", format(nrow(b), big.mark = ","), 100 * b[, mean((h1f_T > 0.5) == (h1f_N > 0.5))])
  msg("    all genes: same major haplotype %.1f%% ; Spearman(h1 freq T, h1 freq N) %.3f", 100 * X[, mean((h1f_T > 0.5) == (h1f_N > 0.5))], X[, cor(h1f_T, h1f_N, method = "s")])
  fl <- list.files(file.path(O, "two"), full.names = TRUE)
  if (length(fl)) { R <- rbindlist(lapply(fl, fread)); pr <- R[, .(t = sum(p_Tfirst < 0.05), n = sum(p_Nfirst < 0.05)), by = pid]
    msg("(2) two-sample (phased), %d patients, %s patient-genes:", uniqueN(R$pid), format(nrow(R), big.mark = ","))
    for (th in c(0.05, 0.01, 0.001)) { st <- R[p_Tfirst < th, .N]; sn <- R[p_Nfirst < th, .N]
      msg("    differential ASE P<%g: tumour-first %s genes ; normal-first %s genes ; ratio %.2f", th, format(st, big.mark = ","), format(sn, big.mark = ","), st / sn) }
    msg("    per patient (P<0.05): tumour-first more %d, normal-first more %d, equal %d ; Wilcoxon P=%.2g", pr[, sum(t > n)], pr[, sum(n > t)], pr[, sum(t == n)], wilcox.test(pr$t, pr$n, paired = TRUE)$p.value) } }
pass("one"); summ(); msg("pass one DONE")
pass("two"); summ(); msg("DONE")
