library(tidyverse)
library(here)
library(ggrastr)
library(ggtext)
library(patchwork)
source(here("code/R/project_defaults_shared.R"))
clr_lab <- c(clr_pheno, mirang = clr_default[2])

# from plot_main_prep_pi
p_pi <- readRDS(here("results/img/R/p_pi.Rds"))

# from plot_main_prep_het
p_het <- readRDS(here("results/img/R/p_het_ind_bp_pheno.Rds"))

# from plot_main_prep_roh
p_froh <- readRDS(here("results/img/R/p_f_rho_callable_ind.Rds"))
p_froh_cum <- readRDS(here("results/img/R/p_cum_f_rho_callable_pheno.Rds"))

# from plot_main_prep_load
p_load_type <- readRDS(here("results/img/R/p_load_by_type_b.Rds")) + theme(axis.ticks = element_line(color ="black"))
p_load_ind <- readRDS(here("results/img/R/p_load_by_ind_b.Rds"))

# from plot_main_prep_load_along_genome
p_load_genome <- readRDS(here("results/img/R/p_load_along_genome.Rds")) +
  guides(color = guide_legend(theme = theme(legend.title.position = "top",
                                            legend.title = element_text(hjust = .5)),
                              override.aes = list(size = 2))) +
  theme(panel.spacing.y = unit(10,"pt"))

set.seed(42)
p_out <- (free(p_pi) + free(p_het) + free(p_froh) + plot_layout(widths = c(.7,1,1))) /
  (p_load_type + theme(axis.title.y = element_text(vjust = .5)) + free(p_load_ind) + plot_layout(widths = c(.8,1))) /
  p_load_genome +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(family = fnt_sel),
        plot.subtitle = element_blank(),
        legend.position = "bottom",
         legend.key.width = unit(11, "pt"),
        legend.key.height = unit(11,"pt"))

# p_out <- cowplot::ggdraw(p_out) + cowplot::draw_label(label = "draft", x = .5, y = .5, hjust = .5, vjust = .5, angle = 30, color = rgb(.6, 0,0,.4), size = 67)

set.seed(42) # fix jitter
ggsave(filename = here("results/img/final/f_whg.pdf"),
       plot = p_out, 
       width = 11,
       height = 10,
       device = cairo_pdf)

set.seed(42) # fix jitter
ggsave(filename = here("results/img/final/f_whg.png"),
       plot = p_out, 
       width = 11,
       height = 10)

p_out_alt <- (free(p_pi) + free(p_het) + free(p_froh) + plot_layout(widths = c(.7,1,1))) /
  (p_load_type + theme(axis.title.y = element_text(vjust = .5)) + free(p_load_ind) + plot_layout(widths = c(.8,1))) /
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(family = fnt_sel),
        plot.subtitle = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(11, "pt"),
        legend.key.height = unit(11,"pt"))

set.seed(42) # fix jitter
ggsave(filename = here("results/img/final/fig4_alt.pdf"),
       plot = p_out_alt, 
       width = 11,
       height = 6.6,
       device = cairo_pdf)

set.seed(42) # fix jitter
ggsave(filename = here("results/img/final/fig4_alt.png"),
       plot = p_out_alt, 
       width = 11,
       height = 6.6)

