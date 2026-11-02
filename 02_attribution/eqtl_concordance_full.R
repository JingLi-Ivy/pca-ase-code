#!/usr/bin/env Rscript
# FAIR eQTL causal test: full CPGEA cis-eQTL anno, best-per-gene (parallel to sqtl_concordance.R).
# Fixes the earlier asymmetry (eQTL used only 66 fig4-candidate eQTLs vs sQTL full 1756).
suppressPackageStartupMessages({ library(data.table) })
PROJ <- "/path/to/data_root/ASE"
EQTL <- "/path/to/data_root/CPGEA/QTLs/CPGEA_eQTL_and_sQTL/CPGEA_cis_tumor_eQTLs_anno.txt"
ASEF <- file.path(PROJ,"all_T_samples_SNPs_with_ASE_annotation.csv")
OUT  <- file.path(PROJ,"figures/eqtl_concordance_full"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
MIN_COUNT <- 10L
lg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), sprintf(...)))

eq <- fread(EQTL)
eq <- eq[FDR < 0.05]
eq[, best := (p.value == min(p.value)), by=gene]
eq <- eq[best==TRUE]
eq <- unique(eq, by="gene")
eq[, chrom := sub("^chr","",as.character(Chr))]
eq[, pos := as.integer(Pos)]
eq[, eqtl_high := ifelse(beta>0, `Minor.Allele`, `Major.Allele`)]
eq <- eq[eqtl_high %in% c("A","C","G","T") & `Minor.Allele`%in%c("A","C","G","T") & `Major.Allele`%in%c("A","C","G","T")]
lg("eQTL best-per-gene (FDR<0.05): %d genes", nrow(eq))

ase <- fread(ASEF)
setnames(ase, old=c("ref.matches","alt.matches"), new=c("ref_n","alt_n"), skip_absent=TRUE)
ase[, chrom := as.character(seqtnames)]; ase[, pos := as.integer(start)]
ase[, n_reads := ref_n + alt_n]
ase <- ase[is.finite(n_reads) & n_reads>=MIN_COUNT & ref_n!=alt_n]
ase[, ase_high := ifelse(ref_n>alt_n, ref, alt)]

m <- merge(ase, eq[,.(chrom,pos,gene,beta,`Minor.Allele`,`Major.Allele`,eqtl_high)],
           by=c("chrom","pos"), allow.cartesian=TRUE)
m <- m[(`Minor.Allele`==ref & `Major.Allele`==alt) | (`Minor.Allele`==alt & `Major.Allele`==ref)]
m[, concordant := (eqtl_high==ase_high)]
lg("Matched eQTL-ASE het obs: %d across %d SNPs, %d genes",
   nrow(m), uniqueN(m[,paste(chrom,pos)]), uniqueN(m$gene))
fwrite(m[,.(chrom,pos,snp_id,gene,beta,eqtl_high,ref,alt,ref_n,alt_n,ase_high,concordant,sample_id)],
       file.path(OUT,"eqtl_full_matched_obs.csv"))

bp <- function(k,n) if(n>0) binom.test(k,n,0.5)$p.value else NA
n_obs<-nrow(m); k_obs<-sum(m$concordant)
snp <- m[, .(ref_n=sum(ref_n), alt_n=sum(alt_n), gene=gene[1], eqtl_high=eqtl_high[1],
             ref=ref[1], alt=alt[1]), by=.(chrom,pos)]
snp[, ase_high:=ifelse(ref_n>alt_n,ref,alt)]; snp<-snp[ref_n!=alt_n]; snp[,concordant:=(eqtl_high==ase_high)]
n_snp<-nrow(snp); k_snp<-sum(snp$concordant)
set.seed(1); N_PERM <- 2000L
perm <- replicate(N_PERM, { fl<-sample(c(TRUE,FALSE),n_snp,replace=TRUE); mean(ifelse(fl,!snp$concordant,snp$concordant)) })
p_perm <- (1+sum(perm >= k_snp/n_snp))/(1+N_PERM)

sink(file.path(OUT,"eqtl_full_concordance_summary.txt"))
cat("eQTL-ASE directional concordance (FULL best-per-gene, parallel to sQTL)\n")
cat(sprintf("eQTL best-per-gene: %d ; matched obs %d (%d SNPs, %d genes)\n\n", nrow(eq), n_obs, n_snp, uniqueN(m$gene)))
cat(sprintf("Per observation : %d/%d = %.1f%%  (binom P=%.3g)\n", k_obs,n_obs,100*k_obs/n_obs, bp(k_obs,n_obs)))
cat(sprintf("Per SNP (pooled): %d/%d = %.1f%%  (binom P=%.3g ; permutation P=%.3g)\n", k_snp,n_snp,100*k_snp/n_snp, bp(k_snp,n_snp), p_perm))
sink()
cat(readLines(file.path(OUT,"eqtl_full_concordance_summary.txt")), sep="\n")

suppressPackageStartupMessages(library(ggplot2))
plt <- data.table(level=c("Per observation","Per SNP"), frac=c(k_obs/n_obs,k_snp/n_snp), n=c(n_obs,n_snp))
p <- ggplot(plt, aes(level,frac))+ geom_col(width=.6, fill="#3B6EA8")+
  geom_hline(yintercept=.5, linetype=2, color="grey40")+
  geom_text(aes(label=sprintf("%.1f%%\n(n=%d)",100*frac,n)), vjust=-0.2, size=3.6)+ ylim(0,1)+
  labs(x=NULL, y="eQTL-ASE directional concordance",
       title="eQTL increasing-allele = RNA high allele",
       subtitle=sprintf("per-SNP %.1f%%; binomial P=%.1e; %d genes",100*k_snp/n_snp,bp(k_snp,n_snp),uniqueN(m$gene)))+
  theme_classic(base_size=11)+theme(plot.subtitle=element_text(size=10))
ggsave(file.path(OUT,"Fig_eqtl_ase_concordance_full.pdf"), p, width=6.2, height=5, device=cairo_pdf)
ggsave(file.path(OUT,"Fig_eqtl_ase_concordance_full.png"), p, width=6.2, height=5, dpi=200)
lg("Done. figure + permutation written.")
