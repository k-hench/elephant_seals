library(tidyverse)
library(here)
library(ggrastr)
library(ggtext)
library(patchwork)
source(here("code/R/project_defaults_shared.R"))
clr_lab <- c(clr_pheno, mirang = clr_default[2])

# from plot_main_prep_pi
p_pi <- readRDS(here("results/img/R/p_pi.Rds"))

# from plot_main_prep_roh
p_froh <- readRDS(here("results/img/R/p_f_rho_callable.Rds"))
p_roh_whg <- readRDS(here("results/img/R/p_rho_whg_callable.Rds"))
p_roh_length <- readRDS(here("results/img/R/p_rho_length_callable.Rds"))

p_froh_cert <- readRDS(here("results/img/R/p_f_rho_certain.Rds"))
p_roh_whg_cert <- readRDS(here("results/img/R/p_rho_whg_certain.Rds"))
p_roh_length_cert <- readRDS(here("results/img/R/p_rho_length_certain.Rds"))

# from plot_main_prep_het
p_het_ind_bp <- readRDS(here("results/img/R/p_het_ind_bp.Rds"))
p_het_ind_snp <- readRDS(here("results/img/R/p_het_ind_snp.Rds"))
p_het_whg_bp <- readRDS(here("results/img/R/p_het_whg_bp.Rds"))
p_het_whg_snp <- readRDS(here("results/img/R/p_het_whg_snp.Rds"))

# from plot_main_prep_load
p_het_load_type <- readRDS(here("results/img/R/p_load_by_type.Rds"))
p_het_load_ind <- readRDS(here("results/img/R/p_load_by_ind.Rds"))

pw1 <- p_pi + theme(axis.line = element_line(linewidth = .35)) + 
  p_froh + theme(axis.line = element_line(linewidth = .35)) +
  p_froh_cert + theme(axis.line = element_line(linewidth = .35))

pw2a <- p_roh_whg + p_roh_length +
  p_roh_whg_cert + p_roh_length_cert +
  plot_layout(widths = c(1,.6))

pw3 <- p_het_ind_bp +
  p_het_whg_bp + 
  p_het_ind_snp +
  p_het_whg_snp +
  plot_layout(nrow = 1)

pw4 <- p_het_load_type + 
  p_het_load_ind +
  plot_layout(widths = c(1,.5))

p_out <- pw1 /
  pw3 /
  pw4 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


p_out2 <- pw2 /
  p_out

ggsave(filename = "results/img/main_all_pannels.pdf",
       plot = p_out2, 
       width = 16,
       height = 22,
       device = cairo_pdf)
