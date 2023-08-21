library(tidyverse)
library(glue)
library(here)

data_all <- read_tsv(here("results", "mutation_load", "snp_eff", "snp_tally","snp_load_pop_details.tsv.gz"))

data_all |>
  filter(!mirang) |>
  mutate(from = pos - 1) |>
  select(chr, from, pos)  |>
  write_tsv(here("results", "mutation_load", "snp_eff", "snp_tally","fixed_in_mirang.bed.gz"), col_names = FALSE)

data_all |>
  filter(!mirleo) |>
  mutate(from = pos - 1) |>
  select(chr, from, pos)  |>
  write_tsv(here("results", "mutation_load", "snp_eff", "snp_tally","fixed_in_mirleo.bed.gz"), col_names = FALSE)
