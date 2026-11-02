suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PUB<-"/path/to/data_root/ASE/figures/pub_figs"
QD<-"/path/to/data_root/CPGEA/QTLs/CPGEA_eQTL_and_sQTL/"
ov<-fread("/path/to/data_root/ASE/code/ASE_SNPs_with_QTL_annotation_exact_match.csv")
ov[,gid:=sub("\\..*","",gene_id)]
pca<-fread("/path/to/data_root/ASE/code/PCa_ASE_genes.csv"); pcag<-unique(sub("\\..*","",pca$gene_id))
ase<-ov[,.(pca=any(gid %in% pcag), gid=gid[1]),by=pos_id]
cat("ASE pos:",nrow(ase)," in PCa-ASE genes:",sum(ase$pca),"\n")
te<-fread(paste0(QD,"CPGEA_cis_tumor_eQTLs_anno.txt")); ne<-fread(paste0(QD,"CPGEA_cis_normal_eQTLs_anno.txt"))
eqpos<-unique(c(paste0(gsub("chr","",te$Chr),":",te$Pos),paste0(gsub("chr","",ne$Chr),":",ne$Pos)))
ts<-fread(paste0(QD,"CPGEA_cis_tumor_sQTLs_anno.txt")); ns<-fread(paste0(QD,"CPGEA_cis_normal_sQTLs_anno.txt"))
sqpos<-unique(c(paste0(gsub("chr","",ts[best_hit==1]$var_chr),":",ts[best_hit==1]$var_from),
                paste0(gsub("chr","",ns[best_hit==1]$var_chr),":",ns[best_hit==1]$var_from)))
ase[,eqtl:=as.integer(pos_id %in% eqpos)][,sqtl:=as.integer(pos_id %in% sqpos)][,qtl:=as.integer(eqtl|sqtl)]
fr<-fread("/path/to/data_root/ASE/derived_inputs/geno_freq_pos.afreq")
setnames(fr,c("#CHROM","ID","REF","ALT","ALT_FREQS","OBS_CT"),c("chrom","id","ref","alt","af","obs"),skip_absent=TRUE)
fr[,maf:=pmin(af,1-af)]; fr<-unique(fr[,.(pos_id=id,maf)],by="pos_id"); ase<-merge(ase,fr,by="pos_id",all.x=TRUE)
fpkm<-fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls")
tc<-grep("^T\\d+_WTS$",names(fpkm),value=TRUE); fpkm[,medF:=apply(.SD,1,median,na.rm=TRUE),.SDcols=tc]
ase<-merge(ase,fpkm[,.(gid=sub("\\..*","",gene_id),medF)],by="gid",all.x=TRUE)
ORf<-function(fg,hit){t<-table(factor(fg,c(FALSE,TRUE)),factor(hit,c(0,1)))
  if(any(dim(t)<2))return(c(OR=NA,lo=NA,hi=NA,p=NA)); f<-fisher.test(t)
  c(OR=unname(f$estimate),lo=f$conf.int[1],hi=f$conf.int[2],p=f$p.value)}
d<-ase[!is.na(maf)]; d[,mbin:=cut(maf,seq(0,0.5,0.05),include.lowest=TRUE)]
res<-list()
for(q in c("eqtl","sqtl")){
  r1<-ORf(d$pca,d[[q]])
  fgd<-d[pca==TRUE]; bgd<-d[pca==FALSE]; need<-fgd[,.N,by=mbin]; setkey(need,mbin); set.seed(1)
  keep<-bgd[,{k=need[.(.BY$mbin),N]; if(is.na(k))k=0L; .SD[sample(.N,min(.N,k))]},by=mbin]
  dm<-rbind(fgd,keep); r2<-ORf(dm$pca,dm[[q]])
  res[[q]]<-rbind(data.table(qtl=q,bg="All tested genes",t(r1)),data.table(qtl=q,bg="MAF-matched",t(r2)))
}
D<-rbindlist(res); D[,qtl:=ifelse(qtl=="eqtl","eQTL","sQTL")]; fwrite(D,file.path(PUB,"S6D_data.csv"))
D[,ylab:=factor(paste(qtl,bg,sep=" | "),levels=rev(paste(rep(c("eQTL","sQTL"),each=2),c("All tested genes","MAF-matched"),sep=" | ")))]
pd<-ggplot(D,aes(OR,ylab,color=qtl))+geom_vline(xintercept=1,linetype=2,color="grey60")+
  geom_pointrange(aes(xmin=lo,xmax=hi),size=.7)+geom_text(aes(label=sprintf("OR=%.3f",OR)),vjust=-.8,size=3,show.legend=FALSE)+
  scale_color_manual(values=c(eQTL="#2166AC",sQTL="#C0392B"),guide="none")+
  labs(x="Odds ratio (PCa-ASE gene SNP enrichment for QTL)",y=NULL)+theme_bw(base_size=12)+theme(panel.grid.minor=element_blank())
ggsave(file.path(PUB,"S6D_ASE_QTL_enrichment_matched.pdf"),pd,width=7,height=3.6)
ggsave(file.path(PUB,"S6D_ASE_QTL_enrichment_matched.png"),pd,width=7,height=3.6,dpi=300)
s<-ase[!is.na(maf)&!is.na(medF)&medF>0]
s[,etert:=cut(medF,quantile(medF,c(0,1/3,2/3,1)),labels=c("Low","Mid","High"),include.lowest=TRUE)]
s[,mafb:=cut(maf,c(0,0.1,0.2,0.3,0.5),labels=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.5"),include.lowest=TRUE)]
strat<-list()
for(et in levels(s$etert))for(mb in levels(s$mafb)){sub<-s[etert==et&mafb==mb]
  if(nrow(sub)>50&&sum(sub$pca)>5){r<-ORf(sub$pca,sub$qtl); strat[[paste(et,mb)]]<-data.table(etert=et,mafb=mb,OR=r["OR"],p=r["p"],n=nrow(sub))}}
S<-rbindlist(strat); S[,etert:=factor(etert,levels=c("Low","Mid","High"))]; fwrite(S,file.path(PUB,"S6F_data.csv"))
pf<-ggplot(S,aes(mafb,etert,fill=OR))+geom_tile(color="white",linewidth=.4)+geom_text(aes(label=sprintf("%.2f",OR)),size=3.6)+
  scale_fill_gradient2(low="#2166AC",mid="white",high="#C0392B",midpoint=1,name="OR")+
  labs(x="Minor allele frequency",y="Gene expression tertile")+theme_minimal(base_size=12)+theme(panel.grid=element_blank())
ggsave(file.path(PUB,"S6F_QTL_enrichment_stratified.pdf"),pf,width=6.2,height=4)
ggsave(file.path(PUB,"S6F_QTL_enrichment_stratified.png"),pf,width=6.2,height=4,dpi=300)
cat("S6D:\n");print(D[,.(qtl,bg,OR=round(OR,2),p)]);cat("S6F strata:",nrow(S),"\nDONE\n")
