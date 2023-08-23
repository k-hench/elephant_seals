library(tidyverse)
library(glue)
library(here)
source(here("code/R/project_defaults.R"))

samples <- read_tsv("data/file_info.tsv") |> 
  select(sample_id, spec, treatment) |> 
  filter(!duplicated(sample_id))

read_all_load <- \(load_type, sample_id){
  data <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}/{sample_id}_{load_type}.bed.gz")))
  data_roh <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}_in_roh/{sample_id}_{load_type}_in_roh.bed.gz")))

  tibble(sample_id = rep(sample_id, 2)) |> 
    left_join(samples) |> 
    mutate(load_type = rep(load_type, 2),
           snp_subset = c("all", "roh"),
           n_snps = c(nrow(data), nrow(data_roh)),
           data = list(data, data_roh))
}

load_types <- c("expressed", "masked", "fixed")
data_load <- expand_grid(sample_id = samples$sample_id,
              load_type = load_types) |> 
  pmap_dfr(read_all_load) |> 
  mutate(load_type = factor(load_type, levels = load_types))

# !! FIXED LOAD NEEDS TO BE CONSTANT WITHIN SPECIES 
# (the snps should be fixed after all)

data_load |> 
  ggplot(aes(x = load_type, y = n_snps, color = treatment)) +
  ggbeeswarm::geom_beeswarm(alpha = .8) +
  # geom_violin() +
  facet_grid(snp_subset ~ spec, scales = "free")
