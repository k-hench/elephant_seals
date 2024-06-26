#!/usr/bin/env Rscript
# Rscript qc_vcfstats.R <ref_tag>
library(tidyverse)
library(here)
library(ggstance)
library(prismatic)
library(glue)
source("code/R/project_defaults.R")

args <- commandArgs(trailingOnly = TRUE)
ref <- as.character(args[1])

files <- c(here(glue("results/genotyping/raw/{ref}_raw_snps_stat.txt")),
           dir(here("results/genotyping/filtered/"),
               pattern = glue("{ref}.*.txt"), full.names = TRUE))

get_vcfstats <- \(file){
  tibble(raw = read_lines(file) |> 
    str_remove("kept ") |> 
    str_replace(" out of [a-z ]*", ",") |> 
    str_replace(" ", ",")) |> 
    separate(raw, into = c("n","n_total", "stat"), convert = TRUE) |> 
    select(-n_total) |> 
    # pivot_wider(names_from = stat,values_from = n) |> 
    mutate(file = str_remove(file, ".*/") |> str_remove("_stat.txt"),
           n = if_else(stat == "Sites", n*1e-6, n),
           stat = if_else(stat == "Sites", "SNPs (n X 10<sup>6</sup>)", "n Individuals"))
  }

p <- files |> 
  map_dfr(get_vcfstats) |> 
  mutate(file = factor(file, levels = c("mirang_raw_snps", "mirang_filtered","mirang_filtered_all",
                                        "mirang_bi-allelic", "mirang_mac1",
                                        "mirang_filtered_mirang", "mirang_filtered_mirleo") |> 
                         rev())) |> 
  ggplot(aes(x = n, y = file)) +
  geom_barh(stat = "identity",
            color = "gray60", fill = clr_alpha("gray85")) +
  geom_text(color = "gray60",
            fill = clr_alpha("gray85"),
            aes(x = ifelse(n < 5, n + .5, .5), 
                label = if_else(n %in% c(20, 40), n, n*1e6)),
            hjust = 0) +
  facet_wrap(stat ~. , scales = "free", nrow = 2) +
  theme_minimal(base_family = fnt_sel) +
  theme(strip.text = ggtext::element_markdown())
  
ggsave(filename = here("results/img/qc/vcf_stats.pdf"),
       plot = p,
       width = 5, height = 5, device = cairo_pdf)
