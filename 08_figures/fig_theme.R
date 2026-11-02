# fig_theme.R — shared publication style for the v9 figures (sourced by the redraw scripts).
# Every panel is saved at its final printed size, so 8 pt here is 8 pt in print:
# Arial (installed in ~/.fonts), all text >= 8 pt including tick labels, no in-panel titles or subtitles
# (information goes to the legend), American spelling ("tumor"), gene symbols in italics (fontface or expression()).
suppressPackageStartupMessages({ library(ggplot2) })
FONT <- "Arial"; BASE <- 8
PAL <- c(tumor = "#B2182B", normal = "#2166AC", grey = "#9E9E9E", dark = "#333333", accent = "#E08214")
theme_pub <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.text = element_text(size = base, colour = "black"), axis.title = element_text(size = base),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          legend.text = element_text(size = base), legend.title = element_text(size = base),
          legend.key.size = unit(3, "mm"), legend.background = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = base),
          plot.title = element_blank(), plot.subtitle = element_blank(),
          plot.margin = margin(2, 2, 2, 2, "mm"))
}
save_panel <- function(p, file, width, height) {
  ggsave(file, p, width = width, height = height, units = "mm", device = cairo_pdf)
  ggsave(sub("\\.pdf$", ".png", file), p, width = width, height = height, units = "mm", dpi = 600)
}
