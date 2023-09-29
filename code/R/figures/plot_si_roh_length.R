library(tidyverse)
library(ggdist)
library(ggtext)
library(here)
library(glue)
library(patchwork)
source(here("code/R/project_defaults_shared.R"))

p_roh <- readRDS(here("results/img/R/p_cum_f_rho_callable_thresholds.Rds"))

p_out <- p_roh +
  plot_layout(guides = "collect") &
  scale_color_manual("Species",
                     values = c(mirang = clr_default[[1]],
                                mirleo = clr_default[[2]]),
                     labels = spec_names,
                     guide = guide_legend(override.aes = list(linewidth = 2,
                                                              alpha = 1)))  &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 11, hjust = 0, vjust = 0,
                                  family = fnt_sel, margin = margin(0,0,10,0)),
          plot.subtitle = element_blank(),
          plot.margin = margin(15,5,5,5),
          legend.position = "bottom",
          legend.text = element_text(family = fnt_sel, 
                                     face = "italic"),
          legend.margin = margin(-5,0,-2,0)) &
  plot_annotation(tag_levels = list(str_c(c("a",
                                            "b",
                                            "c",
                                            "d"),
                                          ") min ROH length: ",
                                          c(1, 10, 100, 1),
                                          c("kb",
                                            "kb",
                                            "kb",
                                            "Mb")))) 

ggsave(plot = p_out,
       filename = here("results/img/final/sf_roh_length.pdf"),
       width = 10, height = 3.5,
       dev = cairo_pdf)

ggsave(plot = p_out,
       filename = here("results/img/final/sf_roh_length.png"),
       width = 10, height = 3.5)

# p_roh &
#   scale_color_manual(values = c(mirang = clr_default[[1]], 
#                               mirleo = clr_default[[2]]),
#                    guide = "none")
