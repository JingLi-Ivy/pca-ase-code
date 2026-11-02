suppressMessages({library(survival);library(survminer);library(ggplot2);library(data.table);library(limma);library(dplyr)})
BB<-"/path/to/data_root/ASE/code/fig4_results"; O<-"/path/to/data_root/ASE/code/tcga_defB"
dir.create(O,showWarnings=FALSE,recursive=TRUE)
load(file.path(BB,"tcga_validation/TCGA_validation_workspace.RData"))
cm<-clinical_matched
gv<-function(p){c<-grep(p,colnames(cm),value=TRUE,ignore.case=TRUE)[1]; if(is.na(c))return(rep(NA,nrow(cm))); cm[[c]]}
GS<-suppressWarnings(as.numeric(if("GS_numeric"%in%colnames(cm)) cm$GS_numeric else gv("^gleason_score$")))
age<-suppressWarnings(as.numeric(gv("age_at_initial_pathologic_diagnosis")))
Tn<-as.numeric(gsub("[^0-9]","",substr(as.character(gv("^pathology_T_stage$")),2,2)))
fu<-as.numeric(cm[["days_to_last_followup"]]); if(all(is.na(fu))) fu<-as.numeric(cm[["patient.days_to_last_followup"]])
bd<-as.numeric(cm[["patient.days_to_first_biochemical_recurrence"]])
st<-tolower(as.character(cm[["patient.biochemical_recurrence"]]))
ev<-ifelse(st%in%c("yes","true","1"),1,ifelse(st%in%c("no","false","0"),0,NA))
tv<-ifelse(ev==1, dplyr::coalesce(bd,fu), fu)
base<-data.frame(id=cm$tcga_full_id,ev=ev,t=tv/30.44,GS=GS,Tnum=Tn,age=age,stringsAsFactors=FALSE)
base<-base[!is.na(base$ev)&!is.na(base$t)&base$t>0&!is.na(base$id),]
cat("n=",nrow(base)," events=",sum(base$ev),"\n")
g<-ifelse(GS>=8,"high",ifelse(GS<8,"low",NA)); k<-!is.na(g)
f<-eBayes(lmFit(expr_matrix_log2[,k],model.matrix(~factor(g[k],levels=c("low","high")))),trend=TRUE,robust=TRUE)
tt<-topTable(f,coef=2,number=Inf,sort.by="none")
de<-data.table(gene_name=rownames(expr_matrix_log2),log2FC=tt$logFC,p_limma=tt$P.Value,fdr_limma=tt$adj.P.Val)
res<-rbindlist(lapply(rownames(expr_matrix_log2),function(gn){
  d<-merge(base,data.frame(id=colnames(expr_matrix_log2),e=as.numeric(expr_matrix_log2[gn,])),by="id")
  d$grp<-factor(ifelse(d$e>=median(d$e,na.rm=TRUE),"High","Low"),levels=c("Low","High")); d$z<-as.numeric(scale(d$e))
  lr<-tryCatch({l<-survdiff(Surv(t,ev)~grp,data=d);1-pchisq(l$chisq,1)},error=function(x)NA_real_)
  s1<-tryCatch(summary(coxph(Surv(t,ev)~grp,data=d)),error=function(x)NULL)
  s2<-tryCatch(summary(coxph(Surv(t,ev)~z+GS+Tnum+age,data=d)),error=function(x)NULL)
  z2<-tryCatch(cox.zph(coxph(Surv(t,ev)~z+GS+Tnum+age,data=d))$table["z","p"],error=function(x)NA_real_)
  data.table(gene_name=gn,logrank_p=lr,
    HR_split=if(is.null(s1))NA_real_ else s1$conf.int[1,"exp(coef)"],
    HR_lo=if(is.null(s1))NA_real_ else s1$conf.int[1,"lower .95"],
    HR_hi=if(is.null(s1))NA_real_ else s1$conf.int[1,"upper .95"],
    HR_perSD_adj=if(is.null(s2))NA_real_ else s2$coef["z","exp(coef)"],
    p_adj=if(is.null(s2))NA_real_ else s2$coef["z","Pr(>|z|)"],zph_p=z2)}))
res[,fdr_logrank:=p.adjust(logrank_p,"BH")][,fdr_adj:=p.adjust(p_adj,"BH")]
cat("log-rank P<.05=",res[logrank_p<0.05,.N]," BH<.10=",res[fdr_logrank<0.10,.N],
    " | adj Cox p<.05=",res[p_adj<0.05,.N]," BH<.10=",res[fdr_adj<0.10,.N],
    " | PH violation=",res[zph_p<0.05,.N],"/",res[is.finite(zph_p),.N],"\n")
s7<-merge(de,res,by="gene_name"); s7[,DE_highconf:=fdr_limma<0.10 & abs(log2FC)>=0.5]
fwrite(s7[order(fdr_limma)],file.path(O,"Supplementary_Table_S7_TCGA_defB.csv"))
six<-c("KIF18B","NEK2","TOP2A","DRD1","GLB1L3","CNN1")
cat("\n=== Six genes ===\n"); print(s7[gene_name%in%six][order(logrank_p),.(gene_name,
  log2FC=round(log2FC,2),fdr_limma=signif(fdr_limma,3),HR=round(HR_split,2),
  CI=paste0(round(HR_lo,2),"-",round(HR_hi,2)),logrank_p=signif(logrank_p,3),
  fdr_logrank=signif(fdr_logrank,3),HR_adj=round(HR_perSD_adj,2),p_adj=signif(p_adj,3))])
for(gn in six){
  d<-merge(base,data.frame(id=colnames(expr_matrix_log2),e=as.numeric(expr_matrix_log2[gn,])),by="id")
  d$Expression<-factor(ifelse(d$e>=median(d$e,na.rm=TRUE),"High","Low"),levels=c("Low","High"))
  s<-summary(coxph(Surv(t,ev)~Expression,data=d)); l<-survdiff(Surv(t,ev)~Expression,data=d); p<-1-pchisq(l$chisq,1)
  gg<-ggsurvplot(survfit(Surv(t,ev)~Expression,data=d),data=d,
    pval=sprintf("log-rank P = %s\nHR = %.2f (95%% CI %.2f-%.2f)",format.pval(p,digits=2,eps=1e-4),
      s$conf.int[1,"exp(coef)"],s$conf.int[1,"lower .95"],s$conf.int[1,"upper .95"]),
    pval.coord=c(0,0.15),palette=c("#2166AC","#C0392B"),risk.table=TRUE,risk.table.height=0.28,
    legend.title=gn,legend.labs=c("Low","High"),xlab="Months since surgery",
    ylab="BCR-free survival",censor.size=2.5,size=0.9,ggtheme=theme_classic(base_size=11))
  gg$plot<-gg$plot+ggtitle(sprintf("%s (TCGA-PRAD; n = %d, events = %d)",gn,nrow(d),sum(d$ev)))
  pdf(file.path(O,sprintf("KM_%s.pdf",gn)),width=5.4,height=6); print(gg); dev.off()
  png(file.path(O,sprintf("KM_%s.png",gn)),width=5.4,height=6,units="in",res=300); print(gg); dev.off()
}
cat("\nDONE ->",O,"\n")
