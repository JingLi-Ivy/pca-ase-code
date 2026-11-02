#!/usr/bin/env Rscript
# Rebuild matched observations from FACETS DNA pileups and RNA allele labels.
# One patient-allele-site record; align REF/ALT explicitly; patient cluster bootstrap.
suppressPackageStartupMessages(library(data.table));setDTthreads(1)
P<-'/path/to/data_root/ASE';O<-file.path(P,'code/v8_round2/revalidation_v10_20261015');F<-'/path/to/data_root/variants/batch_43/facets'
dir.create(file.path(O,'dna_server_only'),showWarnings=FALSE)
cls<-function(t,l)fifelse(t==2&l==1,'balanced_diploid',fifelse(t==2&l==0,'cn_LOH',fifelse(t==1&l==0,'hemizygous_loss',fifelse(t==0,'homozygous_del',fifelse(t>2&l>=1&(t-l)!=l,'allele_specific_gain',fifelse(t>2&(t-l)==l,'balanced_gain',fifelse(t<2,'loss_other','other')))))))
samples<-sub('\\.tsv$','',list.files(file.path(P,'allelic_direction/ase_split'),pattern='^T[0-9]+\\.tsv$'));samples<-samples[file.exists(file.path(F,'results',paste0(samples,'_facets_fit.rds')))];cat('samples',length(samples),'\n')
worker<-function(s){
 r<-fread(file.path(P,'allelic_direction/ase_split',paste0(s,'.tsv')),header=FALSE,col.names=c('chrom','pos','ref','alt','refm','altm'));r[,chrom:=sub('^chr','',as.character(chrom))];nraw<-nrow(r);r<-r[chrom%in%as.character(1:22)&refm+altm>=10&nchar(ref)==1&nchar(alt)==1];nfiltered<-nrow(r);r<-unique(r);nuniq<-nrow(r);bad<-r[,.N,by=.(chrom,pos,ref,alt)][N>1];r<-r[!bad,on=.(chrom,pos,ref,alt)];r[,rna_alt:=altm/(refm+altm)]
 coords<-file.path(O,'dna_server_only',paste0(s,'_coords.tsv'));fwrite(unique(r[,.(chrom,pos)]),coords,sep='\t',col.names=FALSE)
 raw<-file.path(F,'pileup',paste0(s,'.pileup.csv.gz'));stopifnot(file.exists(raw))
 awk<-'NR==FNR {wanted[$1 SUBSEP $2]=1;next} FNR==1 {print;next} {c=$1;sub(/^chr/,"",c);if((c SUBSEP $2) in wanted)print}'
 cmd<-paste('gzip -dc',shQuote(raw),'| awk -F',shQuote('[,\t]'),shQuote(awk),shQuote(coords),'-')
 p<-fread(cmd=cmd);setnames(p,c('Chromosome','Position','Ref','Alt'),c('chrom','pos','dna_ref','dna_alt'));p[,chrom:=sub('^chr','',as.character(chrom))];pn<-nrow(p);p<-unique(p);pbad<-p[,.N,by=.(chrom,pos,dna_ref,dna_alt)][N>1];p<-p[!pbad,on=.(chrom,pos,dna_ref,dna_alt)];p<-p[File2R+File2A>0&File1R+File1A>0];p[,`:=`(vafT_raw=File2A/(File2R+File2A),vafN_raw=File1A/(File1R+File1A))]
 cat(s,'RNA',nrow(r),'pileup',nrow(p),'\n');flush.console()
 f<-readRDS(file.path(F,'results',paste0(s,'_facets_fit.rds')));js<-as.data.table(f$jointseg);cn<-as.data.table(f$out)[,.(seg,tcn,lcn)];rm(f);js[,chrom:=as.character(chrom)];js<-js[het==1&keep==1&chrom%in%as.character(1:22)&paste(chrom,maploc)%chin%paste(r$chrom,r$pos),.(chrom,pos=maploc,vafT,vafN,seg)];js<-merge(js,cn,by='seg');rm(cn);gc()
 j<-merge(p,js,by=c('chrom','pos'),allow.cartesian=TRUE); cat(s,'CN candidates',nrow(j),'T match',sum(abs(j$vafT_raw-j$vafT)<1e-8),'N match',sum(abs(j$vafN_raw-j$vafN)<1e-8),'\n');flush.console();j<-j[abs(vafT_raw-vafT)<1e-8&abs(vafN_raw-vafN)<1e-8];j<-unique(j[,.(chrom,pos,dna_ref,dna_alt,vafT_raw,tcn,lcn)]);jb<-j[,.N,by=.(chrom,pos,dna_ref,dna_alt)][N>1];j<-j[!jb,on=.(chrom,pos,dna_ref,dna_alt)]
 m<-merge(r,j,by=c('chrom','pos'),allow.cartesian=TRUE);candidate<-nrow(m);m[,orientation:=fifelse(ref==dna_ref&alt==dna_alt,'same',fifelse(ref==dna_alt&alt==dna_ref,'swapped','mismatch'))];mis<-m[orientation=='mismatch',.N];sw<-m[orientation=='swapped',.N];m<-m[orientation!='mismatch'];m[,dna_dev:=fifelse(orientation=='same',vafT_raw,1-vafT_raw)-.5];m[,`:=`(rna_dev=rna_alt-.5,cls=cls(tcn,lcn))];m<-unique(m[,.(chrom,pos,ref,alt,cls,dna_dev,rna_dev)]);mb<-m[,.N,by=.(chrom,pos,ref,alt)][N>1];m<-m[!mb,on=.(chrom,pos,ref,alt)];m[,`:=`(conc=sign(dna_dev)==sign(rna_dev),nontie=dna_dev!=0&rna_dev!=0)];imb<-c('cn_LOH','hemizygous_loss','allele_specific_gain','loss_other');m<-m[!is.na(cls) & (!(cls%in%imb)|abs(dna_dev)>=.1)]
 count<-function(z)z[,.(n=.N,k=sum(conc),n_nontie=sum(nontie),k_nontie=sum(conc&nontie)),by=cls]
 pooled<-m[cls%in%imb,.(n=.N,k=sum(conc),n_nontie=sum(nontie),k_nontie=sum(conc&nontie))][,cls:='POOLED_SCNA']
 counts<-rbind(count(m),pooled,fill=TRUE);counts[,sample:=s]
 audit<-data.table(sample=s,rna_raw=nraw,rna_depthfiltered=nfiltered,rna_exact_unique=nuniq,rna_conflict_allele_keys=nrow(bad),dna_position_rows=pn,dna_raw_conflict_allele_keys=nrow(pbad),DNA_CN_conflict_allele_keys=nrow(jb),candidate_crossjoin_rows=candidate,allele_mismatch_crossjoin_rows=mis,swapped_match_rows=sw,retained_unique_allele_observations=nrow(m),final_conflict_allele_keys=nrow(mb))
 list(counts=counts,audit=audit)
}
if('pilot'%in%commandArgs(TRUE))samples<-intersect(samples,'T513')
cores<-min(16L,as.integer(Sys.getenv('NSLOTS','16')))
res<-parallel::mclapply(samples,worker,mc.cores=cores,mc.preschedule=FALSE)
stopifnot(all(vapply(res,is.list,logical(1))))
c<-rbindlist(lapply(res,'[[','counts'));a<-rbindlist(lapply(res,'[[','audit'));fwrite(c,file.path(O,'dna_server_only/patient_class_counts.csv'));fwrite(a,file.path(O,'dna_server_only/patient_checks.csv'))
set.seed(20261015)
out<-rbindlist(lapply(unique(c$cls),function(cl){z<-c[cls==cl];boot<-replicate(5000,{i<-sample.int(nrow(z),replace=TRUE);sum(z$k[i])/sum(z$n[i])});bn<-replicate(5000,{i<-sample.int(nrow(z),replace=TRUE);sum(z$k_nontie[i])/sum(z$n_nontie[i])});data.table(class=cl,n_observations=sum(z$n),patients=nrow(z),concordance=sum(z$k)/sum(z$n),lower=quantile(boot,.025),upper=quantile(boot,.975),n_nontie=sum(z$n_nontie),concordance_nontie=sum(z$k_nontie)/sum(z$n_nontie),nontie_lower=quantile(bn,.025),nontie_upper=quantile(bn,.975),patient_weighted_mean=mean(z$k/z$n))}))
fwrite(out,file.path(O,'results/DNA_RNA_allele_harmonized_summary.csv'));fwrite(a[,lapply(.SD,sum),.SDcols=setdiff(names(a),'sample')],file.path(O,'results/DNA_RNA_harmonization_checks.csv'));print(out);cat('DONE\n')
