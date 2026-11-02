#!/usr/bin/env Rscript
# Delta-ASE pilot: within the same patient and the same heterozygous SNP, is the germline cis effect (oriented by the
# population eQTL increasing allele) weaker, stronger or unchanged in tumor than in matched normal?
# Sites: tumor locally balanced diploid (FACETS tcn2/lcn1, dominant segment), >= MIN reads in both tissues.
# Orientation: CPGEA normal-tissue cis-eQTL, best variant per gene, FDR < 0.05, beta = Minor.Allele effect
# (same convention as eqtl_concordance_full.R). Tumor-tissue eQTL orientation as sensitivity.
# Purity model: s_T = b1*s_N + b2*s_N*purity ; tumor-cell effect relative to normal = b1 + b2 (at purity 1).
suppressPackageStartupMessages(library(data.table)); setDTthreads(8); set.seed(1)
P <- "/path/to/data_root/ASE"; O <- file.path(P, "code/v8_round2/delta_ase"); dir.create(O, showWarnings = FALSE)
MIN <- as.integer(Sys.getenv("MIN_READS", "10"))
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

rdq <- function(f) { e <- fread(f)[FDR < 0.05]; e[, best := p.value == min(p.value), by = gene]; e <- unique(e[best == TRUE], by = "gene")
  e[, `:=`(chrom = sub("^chr", "", as.character(Chr)), pos = as.integer(Pos), eg = sub("\\..*", "", sub("\\|.*", "", gene)))]
  e[, high := ifelse(beta > 0, Minor.Allele, Major.Allele)]
  e[high %in% c("A","C","G","T") & Minor.Allele %in% c("A","C","G","T") & Major.Allele %in% c("A","C","G","T"),
    .(chrom, pos, eg, beta, minor = Minor.Allele, major = Major.Allele, high)] }
Q <- "/path/to/data_root/CPGEA/QTLs/CPGEA_eQTL_and_sQTL/"
eqN <- rdq(paste0(Q, "CPGEA_cis_normal_eQTLs_anno.txt")); eqT <- rdq(paste0(Q, "CPGEA_cis_tumor_eQTLs_anno.txt"))
msg("best-per-gene eQTLs: normal %d, tumor %d", nrow(eqN), nrow(eqT))

cols <- c("seqtnames","start","ref","alt","gene_id","ref.matches","alt.matches","sample_id")
tum <- fread(file.path(P, "all_T_samples_SNPs_with_ASE_annotation.csv"), select = cols)
fs <- Sys.glob(file.path(P, "N*_allele_specific_counts.pipeline.nonzero.txt")); msg("normal per-SNP files: %d", length(fs))
nor <- rbindlist(lapply(fs, function(f) { d <- fread(f, header = FALSE, select = c(1, 2, 3, 4, 6, 7), col.names = c("seqtnames","start","ref","alt","ref.matches","alt.matches"))
  d[, `:=`(ref = gsub("^b'|'$", "", ref), alt = gsub("^b'|'$", "", alt), sample_id = sub("_.*", "", basename(f)))]; d }))
prep <- function(d) d[, .(pid = gsub("\\D", "", sample_id), chrom = as.character(seqtnames), pos = as.integer(start), ref, alt,
                          eg = sub("\\..*", "", gene_id), rn = ref.matches, an = alt.matches)][rn + an >= MIN]
tum <- unique(prep(tum))
nor <- unique(nor[, .(pid = gsub("\\D", "", sample_id), chrom = as.character(seqtnames), pos = as.integer(start), ref, alt, rn = ref.matches, an = alt.matches)][rn + an >= MIN])
m <- merge(tum, nor, by = c("pid","chrom","pos","ref","alt"), suffixes = c("T","N"))
msg("paired SNP-patient sites with >= %d reads in both tissues: %s (patients %d, genes %d)", MIN, format(nrow(m), big.mark = ","), uniqueN(m$pid), uniqueN(m$eg))

fac <- fread("/path/to/data_root/variants/batch_43/facets/annotated_cnv_output.tsv", select = c("sample_id","gene_name","tcn","lcn","num_mark"))
fac[, `:=`(tcn = suppressWarnings(as.numeric(tcn)), lcn = suppressWarnings(as.numeric(lcn)))]
setorder(fac, sample_id, gene_name, -num_mark); fac <- fac[, .SD[1], by = .(sample_id, gene_name)]
sym <- unique(fread(cmd = "awk -F'\\t' '$3==\"gene\"{g=\"\";n=\"\";if(match($9,/gene_id \"[^\"]+\"/))g=substr($9,RSTART+9,RLENGTH-10);if(match($9,/gene_name \"[^\"]+\"/))n=substr($9,RSTART+11,RLENGTH-12);print g\"\\t\"n}' /path/to/data_root/reference/hg38/genes.gtf",
                  header = FALSE, col.names = c("eg","gene_name"))[, eg := sub("\\..*", "", eg)])
fac <- merge(fac, sym, by = "gene_name", allow.cartesian = TRUE)[, pid := gsub("\\D", "", sample_id)]
m <- merge(m, unique(fac[, .(pid, eg, tcn, lcn)]), by = c("pid","eg"), all.x = TRUE)
pu <- fread(file.path(P, "figures/supp_panels/facets_purity_ploidy_per_sample.csv"))[, .(pid = gsub("\\D", "", sample), purity)]
m <- merge(m, pu, by = "pid", all.x = TRUE)
msg("of which tumor balanced diploid: %s ; unevaluable CN: %s", format(m[tcn == 2 & lcn == 1, .N], big.mark = ","), format(m[is.na(tcn) | is.na(lcn), .N], big.mark = ","))
B <- m[tcn == 2 & lcn == 1 & is.finite(purity)]
fwrite(B, file.path(O, "paired_balanced_sites.csv.gz"))

msg("reference fraction (reads pooled): normal %.4f, tumor %.4f", B[, sum(rnN) / sum(rnN + anN)], B[, sum(rnT) / sum(rnT + anT)])

orient <- function(eq, lab) {
  x <- merge(B, eq, by = c("chrom","pos","eg"))
  x <- x[(minor == ref & major == alt) | (minor == alt & major == ref)]
  x[, `:=`(hN = ifelse(high == ref, rnN, anN), hT = ifelse(high == ref, rnT, anT), dN = rnN + anN, dT = rnT + anT)]
  x[, `:=`(sN = hN / dN - 0.5, sT = hT / dT - 0.5)]
  pat <- x[, .(sN = sum(hN) / sum(dN) - 0.5, sT = sum(hT) / sum(dT) - 0.5, n = .N, purity = purity[1]), by = pid]
  boot <- replicate(2000, { p <- pat[sample(.N, replace = TRUE)]; c(p[, weighted.mean(sN, n)], p[, weighted.mean(sT, n)]) })
  mN <- pat[, weighted.mean(sN, n)]; mT <- pat[, weighted.mean(sT, n)]
  r <- boot[2, ] / boot[1, ]
  msg("\n[%s] eQTL-oriented paired sites: %s obs, %d SNPs, %d genes, %d patients", lab, format(nrow(x), big.mark = ","), uniqueN(x[, paste(chrom, pos)]), uniqueN(x$eg), uniqueN(x$pid))
  msg("[%s] mean signed deviation toward eQTL-high allele: normal %+.4f [%+.4f, %+.4f] ; tumor %+.4f [%+.4f, %+.4f]", lab,
      mN, quantile(boot[1, ], .025), quantile(boot[1, ], .975), mT, quantile(boot[2, ], .025), quantile(boot[2, ], .975))
  msg("[%s] tumor/normal ratio %.3f [%.3f, %.3f] ; concordance with eQTL: normal %.1f%%, tumor %.1f%%", lab, mT / mN, quantile(r, .025), quantile(r, .975),
      100 * x[sN != 0, mean(sN > 0)], 100 * x[sT != 0, mean(sT > 0)])
  f <- lm(sT ~ 0 + sN + sN:purity, data = x, weights = pmin(dN, dT))
  cf <- coef(f); bb <- replicate(1000, { ids <- sample(unique(x$pid), replace = TRUE); xx <- rbindlist(lapply(ids, function(i) x[pid == i]))
    sum(coef(lm(sT ~ 0 + sN + sN:purity, data = xx, weights = pmin(dN, dT)))) })
  msg("[%s] purity model: s_T = %.3f*s_N + %.3f*s_N*purity ; tumor-cell relative effect (purity=1) %.3f [%.3f, %.3f]", lab,
      cf[1], cf[2], sum(cf), quantile(bb, .025), quantile(bb, .975))
  x[, pq := cut(purity, quantile(pat$purity, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low","mid","high"))]
  print(x[, .(patients = uniqueN(pid), obs = .N, normal = sum(hN) / sum(dN) - 0.5, tumor = sum(hT) / sum(dT) - 0.5), by = pq][order(pq)][, ratio := tumor / normal][])
  fwrite(x, file.path(O, sprintf("oriented_sites_%s.csv.gz", lab)))
  data.table(orientation = lab, obs = nrow(x), snps = uniqueN(x[, paste(chrom, pos)]), genes = uniqueN(x$eg), patients = uniqueN(x$pid),
             normal = mN, tumor = mT, ratio = mT / mN, ratio_lo = quantile(r, .025), ratio_hi = quantile(r, .975),
             tumor_cell_rel = sum(cf), tcr_lo = quantile(bb, .025), tcr_hi = quantile(bb, .975)) }
R <- rbind(orient(eqN, "normal_eQTL"), orient(eqT, "tumor_eQTL"))
nc <- B[!paste(chrom, pos) %in% c(eqN[, paste(chrom, pos)], eqT[, paste(chrom, pos)])][sample(.N, min(.N, 200000))]
nc[, flip := sample(c(TRUE, FALSE), .N, replace = TRUE)]
nc[, `:=`(sN = ifelse(flip, rnN, anN) / (rnN + anN) - 0.5, sT = ifelse(flip, rnT, anT) / (rnT + anT) - 0.5)]
msg("\nnegative control (non-eQTL sites, random orientation, n=%s): normal %+.4f tumor %+.4f", format(nrow(nc), big.mark = ","), mean(nc$sN), mean(nc$sT))
fwrite(R, file.path(O, "delta_ase_pilot_summary.csv"))
msg("DONE")
