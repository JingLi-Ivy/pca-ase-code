suppressPackageStartupMessages({library(data.table); library(ggplot2); library(org.Hs.eg.db); library(AnnotationDbi)})
setwd("/path/to/data_root/ASE")
OUT<-"figures/pub_figs"; dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
MINN<-5

read_set<-function(files,tissue){
  L<-lapply(files,function(f){d<-fread(f,select=c("gene_id","a")); d[is.finite(a)]})
  x<-rbindlist(L)
  x[,.(mean_ASE=mean(abs(a)), n=.N), by=gene_id][, tissue:=tissue][]
}
tf<-list.files(".",pattern="^T\\d+_ase_results\\.csv$",full.names=TRUE)
nf<-list.files(".",pattern="^N\\d+_ase_results\\.csv$",full.names=TRUE)
cat(sprintf("T files %d  N files %d\n",length(tf),length(nf)))
Tg<-read_set(tf,"Tumor"); Ng<-read_set(nf,"Normal")
m<-merge(Ng[,.(gene_id,Normal=mean_ASE,nN=n)], Tg[,.(gene_id,Tumor=mean_ASE,nT=n)], by="gene_id")
m<-m[nN>=MINN & nT>=MINN]
cat("genes with both tissues n>=",MINN,":",nrow(m),"\n")

ids<-sub("\\..*$","",m$gene_id)
sym<-mapIds(org.Hs.eg.db,keys=ids,column="SYMBOL",keytype="ENSEMBL",multiVals="first")
m[,symbol:=sym[match(sub("\\..*$","",gene_id),names(sym))]]

onco<-fread("/path/to/data_root/reference/oncoKB_cancerGeneList.tsv")
cancer_sym<-unique(onco[["Hugo Symbol"]])
imprinted<-c("DIRAS3","RNU5D-1","TP73","LRRTM1","MYCN","ZDBF2","GPR1","NAP1L5","ERAP2",
  "RHOBTB3","VTRNA2-1","NR3C1","HMMR","PRIM2","ADTRP","FAM50B","PXDC1","LIN28B","AIM1",
  "PLAGL1","HYMAI","SLC22A2","SLC22A3","GRB10","DDC","GLI3","HECW1","HOXA4","PEG10","MAGI2",
  "SGCE","PPP1R9A","TFPI2","CCDC71L","COPG2IT1","MEST","MESTIT1","KLF14","SVOPL","PTPRN2",
  "DLGAP2","PEG13","KCNK9","MYC","ZFAT","ZFAT-AS1","GLIS3","VIM","INPP5F","WT1","WT1-AS",
  "OSBPL5","KCNQ1DN","PHLDA2","IGF2","SLC22A18AS","KCNQ1","ASCL2","CDKN1C","KCNQ1OT1",
  "IGF2AS","H19","SLC22A18","INS","ANO1","ZC3H12C","NTM","ST8SIA1","RBP5","SMOC1","MEG3",
  "DIO3","DLK1","MEG8","SNORD114-1","SNORD113-1","DIO3OS","RTL1","MAGEL2","MKRN3","UBE3A",
  "SNORD107","SNORD116","SNRPN","SNORD109B","SNORD109A","ATP10A","SNORD108","PWAR6","NPAP1",
  "NDN","SNORD64","SNURF","RASGRF1","IRAIN","PRR25","ZNF597","CMTM1","ZFP90","ZNF396",
  "TCEB3C","PARD6G-AS1","ADNP2","PARD6G","DNMT1","AXL","ZIM2","PEG3","PEG3-AS1","MIMT1",
  "NNAT","BLCAP","MCTS2","GDAP1L1","SGK2","GNAS","L3MBTL1","GNAS-AS1","SIM2","DSCAM",
  "DGCR6","DGCR6L")

mk<-function(d,lab) d[, .(set=lab, Normal, Tumor, gene_id, symbol)]
S_all<-mk(m,"All genes")
S_can<-mk(m[symbol %in% cancer_sym],"Cancer genes")
S_imp<-mk(m[symbol %in% imprinted],"Imprinted controls")
cat(sprintf("set sizes: all %d | cancer %d | imprinted %d\n",nrow(S_all),nrow(S_can),nrow(S_imp)))
cat("imprinted present:",paste(sort(S_imp$symbol),collapse=","),"\n")

allsets<-rbindlist(list(S_all,S_can,S_imp))
fwrite(allsets,file.path(OUT,"Fig1D_threeset_burden_data.csv"))

long<-melt(allsets,id.vars=c("set","gene_id","symbol"),measure.vars=c("Normal","Tumor"),
           variable.name="tissue",value.name="ase")
long[,set:=factor(set,levels=c("All genes","Cancer genes","Imprinted controls"))]
long[,tissue:=factor(tissue,levels=c("Normal","Tumor"))]

stat<-allsets[,{w<-wilcox.test(Tumor,Normal,paired=TRUE);
  .(p=w$p.value, n=.N, ytop=max(c(Normal,Tumor)))}, by=set]
stat[,set:=factor(set,levels=c("All genes","Cancer genes","Imprinted controls"))]
stat[,lab:=ifelse(p<2.2e-16,"P < 2.2e-16",ifelse(p<0.001,sprintf("P = %.1e",p),sprintf("P = %.3f",p)))]
stat[,xnum:=as.integer(set)]
print(stat[,.(set,n,p,lab)])

ymax<-max(long$ase)
TH<-theme_classic(base_size=12)+theme(plot.title=element_blank(),
     legend.position="top", legend.title=element_blank(),
     legend.background=element_blank(), plot.margin=margin(7,14,7,7))

p<-ggplot(long,aes(set,ase,fill=tissue))+
  geom_boxplot(width=.62,outlier.size=.3,position=position_dodge(.7))+
  scale_fill_manual(values=c(Normal="#4DAF4A",Tumor="#E41A1C"))+
  geom_segment(data=stat,aes(x=xnum-.18,xend=xnum+.18,y=ymax*1.04,yend=ymax*1.04),inherit.aes=FALSE,linewidth=.4)+
  geom_text(data=stat,aes(x=xnum,y=ymax*1.10,label=lab),inherit.aes=FALSE,size=3.3)+
  geom_text(data=stat,aes(x=xnum,y=-0.02,label=paste0("n=",n)),inherit.aes=FALSE,size=3,color="grey35")+
  scale_y_continuous(expand=expansion(mult=c(.04,.10)))+
  labs(x=NULL,y="Mean ASE magnitude per gene")+TH

ggsave(file.path(OUT,"Fig1D_tumor_vs_normal_ASE_burden.pdf"),p,width=5.6,height=4.4)
ggsave(file.path(OUT,"Fig1D_tumor_vs_normal_ASE_burden.png"),p,width=5.6,height=4.4,dpi=300)
cat("DONE\n")
