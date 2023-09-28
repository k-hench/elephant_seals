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

set.seed(42)
p_out <- (p_pi + p_het + p_froh) /
  (p_load_type + p_load_ind) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(family = fnt_sel),
        plot.subtitle = element_blank())

ggsave(filename = here("results/img/final/f_whg.pdf"),
       plot = p_out, 
       width = 10,
       height = 6.66,
       device = cairo_pdf)

ggsave(filename = here("results/img/final/f_whg.png"),
       plot = p_out, 
       width = 10,
       height = 6)

p_out_alt <- (p_pi + p_het + p_froh_cum) /
  (p_load_type + p_load_ind) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = fnt_sel),
        plot.subtitle = element_blank())


ggsave(filename = here("results/img/main_draft_alt.pdf"),
       plot = p_out_alt, 
       width = 10,
       height = 6.66,
       device = cairo_pdf)

