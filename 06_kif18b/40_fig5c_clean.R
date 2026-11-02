suppressPackageStartupMessages({library(data.table); library(ggplot2)})
OUT <- "/path/to/data_root/ASE/figures/pub_figs"
d <- fread("/path/to/data_root/ASE/figures/FOXM1_analysis/FOXM1_KIF18B_ASE_frequency.csv")
d[, FOXM1_group := factor(FOXM1_group, levels=c("Low","High"))]
d[, freq_percent := 100*freq_KIF18B_ASE]
m <- matrix(c(d[FOXM1_group=="Low",n_KIF18B_ASE], d[FOXM1_group=="Low",n_total-n_KIF18B_ASE],
              d[FOXM1_group=="High",n_KIF18B_ASE], d[FOXM1_group=="High",n_total-n_KIF18B_ASE]), nrow=2)
fp <- fisher.test(m)$p.value
p <- ggplot(d, aes(FOXM1_group, freq_percent, fill=FOXM1_group)) +
  geom_col(width=0.6) +
  geom_text(aes(label=sprintf("%.1f%%\n(%d/%d)", freq_percent, n_KIF18B_ASE, n_total)), vjust=-0.3, size=3.4) +
  scale_fill_manual(values=c(Low="#00AFBB", High="#FC4E07"), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(title="KIF18B ASE frequency by FOXM1 expression",
       caption=sprintf("Fisher's exact test P = %.3g", fp),
       x="FOXM1 expression group", y="KIF18B ASE frequency, percent") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        plot.caption=element_text(hjust=0.5, size=9, color="grey30"))
ggsave(file.path(OUT,"Fig5C_FOXM1_high_low.pdf"), p, width=4.6, height=4.2)
ggsave(file.path(OUT,"Fig5C_FOXM1_high_low.png"), p, width=4.6, height=4.2, dpi=300)
cat(sprintf("DONE Fisher P=%.4g\n", fp))
