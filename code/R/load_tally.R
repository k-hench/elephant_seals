library(tidyverse)
library(prismatic)
library(patchwork)
library(glue)
library(here)

import_genomes <- function(spec){
  read_tsv(here("data", "genomes", "filtered", glue("{spec}_filt.fa.gz.fai")),
           col_names = c("chr", "length", "offset", "linebases", "linewidth")) |> 
    mutate(spec = spec,
           idx = row_number(),
           end_pos = cumsum(length),
           start_pos = lag(end_pos, default = 0),
           eo = idx %% 2)
}

genome <- import_genomes("mirang") |> 
  select(`#CHROM` = chr, gstart = start_pos)

import_tsv <- \(base){
  read_tsv(here("results", "mutation_load", "snp_eff", "snp_tally", glue("{base}.tsv.gz"))) |> 
    left_join(genome) |> 
    mutate(gpos = gstart + POS,
           new_col = TRUE) |> 
    set_names(nm = c("chr", "pos", "gstart", "gpos", base)) |> 
    select(-c(chr:gstart))
}

# create master-table
# (classifying all SNPs with respect to variance within populations
#  and load status)
data_all <- import_tsv("all") |> 
  left_join(import_tsv("mirang"))|> 
  left_join(import_tsv("mirleo"))|> 
  left_join(import_tsv("load")) |> 
  mutate(across(mirang:load, \(x){replace_na(x, FALSE)}),
         both = mirang & mirleo,
         variable_in = case_when(
           both ~ "both",
           mirang & !mirleo ~ "ang_only",
           !mirang & mirleo ~ "leo_only",
           !mirang & !mirleo ~ "neither"))

# summaise master-table and export
data_all |> 
  group_by(variable_in, load) |> 
  count() |> 
  mutate(load = c("no_load", "load")[load + 1]) |> 
  pivot_wider(names_from = load, values_from = n) |> 
  mutate(permil = load / (load + no_load) * 1000,
         variable_in = factor(variable_in,
                              levels = c("both",
                                         "ang_only",
                                         "leo_only",
                                         "neither",
                                         "Total"))) |> 
  arrange(variable_in) |> 
  as_tibble() |> 
  write_tsv(here("results", "mutation_load", "snp_eff", "snp_tally","n_snp_load_in_pop.tsv"))

# export master-table
data_all|>
  write_tsv(here("results", "mutation_load", "snp_eff", "snp_tally","snp_load_pop_details.tsv.gz"))
