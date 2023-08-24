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
  if(load_type == "masked"){
    data_anc <- data
  } else {
    data_anc <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}_anc/{sample_id}_{load_type}_anc.bed.gz")))
  }

  tibble(sample_id = rep(sample_id, 3)) |> 
    left_join(samples) |> 
    mutate(load_type = rep(load_type, 3),
           snp_subset = c("all", "anc","roh"),
           n_snps = c(nrow(data), nrow(data_anc), nrow(data_roh)),
           data = list(data, data_anc, data_roh))
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
  ggbeeswarm::geom_beeswarm(alpha = .8)+#, aes(color = sample_id == "160488")) +
  # geom_violin() +
  facet_grid(snp_subset ~ spec, scales = "free") +
  # ggrepel::geom_text_repel(mapping = aes(label = sample_id)) +
  theme_bw()


ggsave("~/Downloads/load_messup.pdf", width = 8, height = 8, device = cairo_pdf)
