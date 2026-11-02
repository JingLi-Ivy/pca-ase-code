#!/usr/bin/env Rscript
# 26e: 26b with the correct depth for normal tissue. 26b (and tn_hethom/13/16/26c/26d) merged both tissues with
# ase_gene_sample_agg.csv.gz, which holds TUMOUR samples only, so normal |a| was filtered (dep>=20) and depth-residualised
# on the tumour depth of the same patient. Here normal uses cistrans_full/ase_gene_sample_agg_N.csv.gz. Nothing else changed.
# 26b (review follow-up): slope-interaction version of 13 on PAIRED patients only, with a strict patient permutation of the
# tumour-normal difference. Statistic per gene x TF: y = rank-inverse-normal of the |a| residual pooled over het + hom (within tissue);
# z = (b_het - b_hom) / sqrt(se_het^2 + se_hom^2) from OLS slopes on the TF z-score; unequal noise in y does not bias b.
# Original header of 13:
# v8_round2 (a): is the tumor-normal difference in the heterozygote-specific TF-|a| association explained by tissue composition?
# Same data, variant selection (seed 42, identical code path) and statistic as tn_hethom.R; adds covariate residualization
# with definitions applied identically in tumor and normal:
#   epithelial marker score = mean z of log2 FPKM of EPCAM, KRT8, KRT18, CDH1, CLDN4, KRT19 (within tissue)
#   proliferation score     = mean z of MKI67, PCNA, CCNB1, CCNB2, CDK1, CCNA2, AURKA (as in foxm1_interaction_robustness.R)
# Optional: if env DECON_T / DECON_N point to BayesPrism cell_fraction tables with an Epithelial_cells column, a third model uses them.
suppressPackageStartupMessages({library(data.table); library(parallel)}); setDTthreads(1)
NREP<-as.integer(Sys.getenv("NREP","100")); NC<-as.integer(Sys.getenv("NC","16"))
msg<-function(...){cat(sprintf(...),"\n"); flush.console()}
OUT<-"/path/to/data_root/ASE/code/v8_round2/cistrans_slope_26e"; dir.create(OUT,showWarnings=FALSE)
AD<-"/path/to/data_root/ASE/ADASTRA6/ADASTRA.v.6.1.Mabel/TF"
PP<-"/path/to/data_root/CPGEA/CPGEA_genotype/combined.134.qc"
WIN<-10000; set.seed(42)

bim<-fread(paste0(PP,".bim"),col.names=c("chr","snp","cm","pos","a1","a2"))[,idx:=.I][,chr:=paste0("chr",chr)]
fam<-fread(paste0(PP,".fam"),col.names=c("fid","iid","f","m","sex","pheno"))
sids<-as.character(sub("^N","",fam$fid)); N<-nrow(fam); bps<-ceiling(N/4)
CON<-file(paste0(PP,".bed"),"rb"); readBin(CON,"raw",n=3)
LUT<-matrix(NA_integer_,256L,4L)
for(bv in 0:255) for(j in 0:3){cd<-bitwAnd(bitwShiftR(bv,2L*j),3L)
  LUT[bv+1L,j+1L]<-if(cd==0L)2L else if(cd==1L)NA_integer_ else if(cd==2L)1L else 0L}
readsnp<-function(idx){seek(CON,3+(idx-1)*bps); rb<-readBin(CON,"raw",n=bps)
  as.vector(t(LUT[as.integer(rb)+1L,,drop=FALSE]))[seq_len(N)]}

ann<-unique(fread("/path/to/data_root/ASE/all_T_samples_SNPs_with_ASE_annotation.csv",select=c("seqtnames","start","gene_id")))
ann[,chr:=paste0("chr",as.character(seqtnames))]
gpos<-ann[,.(gchr=chr[1],gmin=min(start),gmax=max(start)),by=gene_id]
fpkm<-fread("/path/to/data_root/CPGEA/CPGEA_web/RNA-seq/FPKM/filter_genes.FPKM.xls")
tcc<-grep("^T\\d+_WTS$",names(fpkm),value=TRUE); ncc<-grep("^N\\d+_WTS$",names(fpkm),value=TRUE)
esT<-as.character(gsub("\\D","",tcc)); esN<-as.character(gsub("\\D","",ncc))
PAIR<-intersect(esT,esN)
ag<-fread("/path/to/data_root/ASE/code/ase_gene_sample_agg.csv.gz")[,sid:=as.character(gsub("\\D","",sample_id))]
agN<-fread("/path/to/data_root/ASE/code/cistrans_full/ase_gene_sample_agg_N.csv.gz")[,sid:=as.character(sid)]

score<-function(genes,cc,es){ m<-log2(as.matrix(fpkm[gene_name %in% genes, ..cc])+1)
  data.table(sid=es, v=as.numeric(colMeans(t(scale(t(m))),na.rm=TRUE))) }
EPI<-c("EPCAM","KRT8","KRT18","CDH1","CLDN4","KRT19"); PRO<-c("MKI67","PCNA","CCNB1","CCNB2","CDK1","CCNA2","AURKA")
covT<-merge(score(EPI,tcc,esT)[,.(sid,epi=v)], score(PRO,tcc,esT)[,.(sid,pro=v)], by="sid")
covN<-merge(score(EPI,ncc,esN)[,.(sid,epi=v)], score(PRO,ncc,esN)[,.(sid,pro=v)], by="sid")
rdec<-function(f){ if(!nzchar(f)||!file.exists(f)) return(NULL); d<-fread(f); d[,.(sid=gsub("\\D","",get(names(d)[1])), bp_epi=Epithelial_cells)] }
dT<-rdec(Sys.getenv("DECON_T")); dN<-rdec(Sys.getenv("DECON_N"))
if(!is.null(dT)) covT<-merge(covT,dT,by="sid",all.x=TRUE); if(!is.null(dN)) covN<-merge(covN,dN,by="sid",all.x=TRUE)
HAS_BP<-!is.null(dT)&&!is.null(dN)
if(HAS_BP){ msg("BayesPrism epithelial fraction: tumor %d, normal %d samples", sum(is.finite(covT$bp_epi)), sum(is.finite(covN$bp_epi)))
  msg("  Spearman(marker score, BayesPrism epi): tumor %.2f, normal %.2f", cor(covT$epi,covT$bp_epi,method="s",use="c"), cor(covN$epi,covN$bp_epi,method="s",use="c")) }

loadase<-function(pre,es,cov){
  fs<-list.files("/path/to/data_root/ASE",pattern=sprintf("^%s\\d+_ase_results\\.csv$",pre),full.names=TRUE)
  a<-rbindlist(lapply(fs,function(f){d<-fread(f)[,.(gene_id,a)];d[,sid:=as.character(gsub("\\D","",basename(f)))];d}))
  d<-merge(a,(if(pre=="N") agN else ag)[,.(gene_id,sid,dep)],by=c("gene_id","sid"))[dep>=20 & sid %in% es]
  d[,`:=`(ai=abs(a), ld=log10(dep))]; merge(d,cov,by="sid",all.x=TRUE)}
aT<-loadase("T",esT,covT)[sid %in% PAIR]; aN<-loadase("N",esN,covN)[sid %in% PAIR]
msg("paired patients with ASE in both tissues: tumour %d, normal %d", uniqueN(aT$sid), uniqueN(aN$sid))
gT<-aT[,.N,by=gene_id][N>=80]$gene_id; gN<-aN[,.N,by=gene_id][N>=80]$gene_id
gok<-intersect(gT,gN); aT<-aT[gene_id %in% gok]; aN<-aN[gene_id %in% gok]
msg("genes measurable in both tissues: %s", format(length(gok),big.mark=","))
res<-function(d,fml){ d[,{ok<-complete.cases(.SD); v<-rep(NA_real_,.N)
  if(sum(ok)>20) v[ok]<-as.numeric(residuals(lm(fml,data=.SD[ok]))); v}, by=gene_id, .SDcols=all.vars(fml)]$V1 }
MODELS<-list(r0=ai~ld, r_epi=ai~ld+epi, r_pro=ai~ld+pro, r_all=ai~ld+epi+pro)
if(HAS_BP) MODELS$r_bp<-ai~ld+bp_epi
for(nm in names(MODELS)){ setorder(aT,gene_id); setorder(aN,gene_id)
  aT[,(nm):=res(aT,MODELS[[nm]])]; aN[,(nm):=res(aN,MODELS[[nm]])] }
setkey(aT,gene_id); setkey(aN,gene_id)

gvv<-function(s,cc){r<-fpkm[gene_name==s]; if(!nrow(r))return(NULL)
  if(nrow(r)>1){m<-apply(as.matrix(r[,..cc]),1,median); r<-r[which.max(m)]}; as.numeric(r[1,..cc])}
map<-fread("/path/to/data_root/ASE/code/adastra_tf_map.csv")
map<-map[file.exists(file.path(AD,paste0(mn,"_HUMAN.tsv")))]
map[,`:=`(mT=sapply(sym,function(s){v<-gvv(s,tcc); if(is.null(v))NA_real_ else median(v,na.rm=TRUE)}),
          mN=sapply(sym,function(s){v<-gvv(s,ncc); if(is.null(v))NA_real_ else median(v,na.rm=TRUE)}))]
tfl<-map[is.finite(mT)&is.finite(mN)&mT>=1&mN>=1]
TZT<-do.call(cbind,lapply(tfl$sym,function(s) as.numeric(scale(log2(gvv(s,tcc)+0.01)))))
TZN<-do.call(cbind,lapply(tfl$sym,function(s) as.numeric(scale(log2(gvv(s,ncc)+0.01)))))
colnames(TZT)<-colnames(TZN)<-tfl$sym; rownames(TZT)<-esT; rownames(TZN)<-esN
int<-function(y) qnorm((rank(y)-0.5)/length(y))
slz<-function(y1,M1,y0,M0){
  f<-function(y,M){n<-length(y); Mc<-sweep(M,2,colMeans(M)); yc<-y-mean(y); sxx<-colSums(Mc^2)
    b<-as.numeric(crossprod(Mc,yc))/sxx; rss<-sum(yc^2)-b^2*sxx; list(b=b,v=pmax(rss,0)/(n-2)/sxx)}
  a<-f(y1,M1); c<-f(y0,M0); z<-(a$b-c$b)/sqrt(a$v+c$v); z[!is.finite(z)]<-0; z}
fz<-function(r1,M1,r0,M0,n1,n0){
  c1<-suppressWarnings(as.numeric(cor(r1,M1,method="spearman"))); c0<-suppressWarnings(as.numeric(cor(r0,M0,method="spearman")))
  c1[!is.finite(c1)]<-0; c0[!is.finite(c0)]<-0
  (atanh(pmin(pmax(c1,-.99),.99))-atanh(pmin(pmax(c0,-.99),.99)))/sqrt(1/(n1-3)+1/(n0-3))}

both<-intersect(esT,esN); setkey(bim,chr); sel<-list()
for(g in gok){
  gp<-gpos[gene_id==g]; if(!nrow(gp)) next
  cd<-bim[.(gp$gchr)][pos>=gp$gmin-WIN & pos<=gp$gmax+WIN]; if(!nrow(cd)) next
  if(nrow(cd)>40) cd<-cd[sample(.N,40)]
  best<-NULL; bv<--1L
  for(k in seq_len(nrow(cd))){H<-setNames(readsnp(cd$idx[k]),sids)[both]
    v<-min(sum(H==1,na.rm=TRUE), sum(H %in% c(0,2),na.rm=TRUE)); if(v>bv){bv<-v; best<-cd[k]}}
  if(is.null(best)||bv<25) next
  sel[[length(sel)+1]]<-list(gene_id=g,idx=best$idx,snp=best$snp,minhh=bv)}
S<-rbindlist(sel); msg("selected genes = %s", format(nrow(S),big.mark=","))
GT<-sapply(S$idx,readsnp); rownames(GT)<-sids

GDT<-lapply(S$gene_id,function(g) aT[.(g)]); GDN<-lapply(S$gene_id,function(g) aN[.(g)])
scan<-function(GDx,TZ,col,G=GT){ out<-vector("list",nrow(S))
  for(i in seq_len(nrow(S))){ z<-GDx[[i]]; y<-z[[col]]; hh<-G[z$sid,i]
    kp<-!is.na(hh)&is.finite(y); i1<-which(kp&hh==1L); i0<-which(kp&hh!=1L); n1<-length(i1); n0<-length(i0)
    if(n1<25||n0<25) next
    yy<-rep(NA_real_,length(y)); yy[kp]<-int(y[kp])
    M<-TZ[z$sid,,drop=FALSE]; out[[i]]<-data.table(g=S$gene_id[i], z=slz(yy[i1],M[i1,,drop=FALSE],yy[i0],M[i0,,drop=FALSE])) }
  rbindlist(out)[,k:=seq_len(.N),by=g] }
gboot<-function(v){ b<-replicate(5000,mean(sample(v,replace=TRUE))); c(mean(v),quantile(b,c(.025,.975)),mean(b<=0)) }
R<-list()
for(nm in names(MODELS)){
  zt<-scan(GDT,TZT,nm); zn<-scan(GDN,TZN,nm)
  D<-merge(zt,zn,by=c("g","k"),suffixes=c("T","N")); D[,d:=zT-zN]
  gm<-D[,.(t=mean(zT),n=mean(zN),d=mean(d)),by=g]
  bt<-gboot(gm$t); bn<-gboot(gm$n); bd<-gboot(gm$d)
  R[[nm]]<-data.table(model=deparse(MODELS[[nm]]), genes=nrow(gm), tests=nrow(D),
    tumor=bt[1],tumor_lo=bt[2],tumor_hi=bt[3], normal=bn[1],normal_lo=bn[2],normal_hi=bn[3], delta=bd[1],delta_lo=bd[2],delta_hi=bd[3],delta_P=bd[4])
  msg("%-22s genes %s | tumor %+.4f [%+.4f,%+.4f] | normal %+.4f [%+.4f,%+.4f] | delta %+.4f [%+.4f,%+.4f]",
      deparse(MODELS[[nm]]), format(nrow(gm),big.mark=","), bt[1],bt[2],bt[3], bn[1],bn[2],bn[3], bd[1],bd[2],bd[3]) }
fwrite(rbindlist(R), file.path(OUT,"tn_slope_by_covariate_model.csv"))
PM<-names(MODELS)[names(MODELS) %in% c("r0","r_bp")]
msg("\nstrict permutation x%d on %d cores (models: %s)", NREP, NC, paste(PM, collapse=", "))
obs<-rbindlist(R)[, .(model=names(MODELS), tumor, normal, delta)]
null<-rbindlist(mclapply(seq_len(NREP), function(b){ set.seed(5000+b); p<-sample.int(N); G<-GT[p,,drop=FALSE]; rownames(G)<-sids
  rbindlist(lapply(PM, function(nm){ zt<-scan(GDT,TZT,nm,G); zn<-scan(GDN,TZN,nm,G)
    D<-merge(zt,zn,by=c("g","k"),suffixes=c("T","N")); gm<-D[,.(t=mean(zT),n=mean(zN),d=mean(zT-zN)),by=g]
    data.table(rep=b, model=nm, tumor=mean(gm$t), normal=mean(gm$n), delta=mean(gm$d)) })) }, mc.cores=NC))
fwrite(null, file.path(OUT,"tn_slope_strict_null.csv"))
for(nm in PM){ o<-obs[model==nm]; nl<-null[model==nm]
  for(v in c("tumor","normal","delta")){ x<-nl[[v]]; ov<-o[[v]]
    msg("%-5s %-6s obs %+.4f | null %+.4f +- %.4f | %.1f SD | two-sided P=%.4f", nm, v, ov, mean(x), sd(x), (ov-mean(x))/sd(x), (1+sum(abs(x-mean(x))>=abs(ov-mean(x))))/(length(x)+1)) } }
msg("DONE")
