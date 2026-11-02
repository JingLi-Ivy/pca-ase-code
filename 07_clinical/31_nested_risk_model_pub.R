#!/usr/bin/env Rscript
# Publication-quality redraw of the nested prognostic model figure: Clinical / Clinical+Expression / Clinical+Expression+ASE.
# Three models in three colors, shared legend on the right, no bracketed titles or subtitles, A/B/C panel tags, PDF + 300-DPI PNG.
# Reads existing results in nested_risk_model_results/ (the bootstrap is not rerun).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(survival)})
RES <- "/path/to/data_root/ASE/nested_risk_model_results"
OUT <- "/path/to/data_root/ASE/figures/pub_figs"

LV  <- c("Clinical","Clinical + Expression","Clinical + Expression + ASE")
COL <- c("Clinical"="#7F7F7F","Clinical + Expression"="#E69F00","Clinical + Expression + ASE"="#C0392B")
relab <- function(x){ x<-as.character(x)
  fifelse(grepl("^M1",x),LV[1], fifelse(grepl("^M2",x),LV[2], LV[3])) }

TH <- theme_classic(base_size=12)+theme(plot.title=element_blank(),
        legend.title=element_blank(), legend.position="right")

ci <- fread(file.path(RES,"Cindex.csv"))
ci[, model:=factor(relab(model),levels=LV)]
pA <- ggplot(ci, aes(model, c_apparent, fill=model))+
  geom_col(width=.7)+
  geom_hline(yintercept=.5, linetype=3, color="grey40")+
  scale_fill_manual(values=COL, drop=FALSE)+
  coord_cartesian(ylim=c(0,.82))+
  labs(x=NULL, y="C-index")+TH+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

ta <- fread(file.path(RES,"timeAUC.csv")); ta[, model:=factor(relab(model),levels=LV)]
pB <- ggplot(ta, aes(time, AUC, color=model))+
  geom_line(linewidth=.8)+geom_point(size=2)+
  geom_hline(yintercept=.5, linetype=3, color="grey40")+
  scale_color_manual(values=COL, drop=FALSE)+
  coord_cartesian(ylim=c(.45,.85))+
  labs(x="Months", y="Time-dependent AUC")+TH

nd0 <- fread(file.path(RES,"nested_model_data.csv")); S0 <- Surv(nd0$surv_time, nd0$bcr_event)
fits <- list(coxph(S0~age+psa+gleason+path_t, data=nd0),
             coxph(S0~age+psa+gleason+path_t+expr_burden_z, data=nd0),
             coxph(S0~age+psa+gleason+path_t+expr_burden_z+ase_burden_z, data=nd0))
U <- 26L; K <- 4L
calib_one <- function(fit, lab){
  su <- summary(survfit(fit, newdata=nd0), times=U, extend=TRUE)$surv
  risk <- 1 - as.numeric(su)
  g <- cut(risk, breaks=quantile(risk, seq(0,1,length.out=K+1)), include.lowest=TRUE, labels=FALSE)
  rbindlist(lapply(sort(unique(g)), function(k){
    s <- summary(survfit(Surv(surv_time,bcr_event)~1, data=nd0[g==k]), times=U, extend=TRUE)
    data.table(model=lab, pred=mean(risk[g==k]), obs=1-s$surv, obs_lo=1-s$upper, obs_hi=1-s$lower)
  }))
}
ca <- rbindlist(Map(calib_one, fits, LV)); ca[, model:=factor(model, levels=LV)]
pdg <- position_dodge(.006)
pC <- ggplot(ca, aes(pred, obs, color=model))+
  geom_abline(slope=1, intercept=0, linetype=2, color="grey50")+
  geom_line(linewidth=.5, position=pdg)+
  geom_errorbar(aes(ymin=obs_lo, ymax=obs_hi), width=.012, position=pdg)+
  geom_point(size=2.2, position=pdg)+
  scale_color_manual(values=COL, drop=FALSE)+
  coord_cartesian(xlim=c(0,.6), ylim=c(0,.66))+
  labs(x="Predicted risk", y="Observed BCR risk")+TH

save2 <- function(p,name,w,h){ ggsave(file.path(OUT,paste0(name,".pdf")),p,width=w,height=h)
  ggsave(file.path(OUT,paste0(name,".png")),p,width=w,height=h,dpi=300) }
save2(pA, "Fig_nested_Cindex",      5.8, 4.0)
save2(pB, "Fig_nested_timeAUC",     6.0, 4.0)
save2(pC, "Fig_nested_calibration", 6.0, 4.0)
cat("DONE\n")
