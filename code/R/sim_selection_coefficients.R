library(tidyverse)
library(glue)
library(here)
library(prismatic)
source(here("code/R/project_defaults.R"))

read_s <- \(idx){
  read_tsv(here(glue("results/simulations/selection_coefficient/sim_{idx}_s.tsv.gz"))) |> 
    mutate(sim_idx = idx)
}

data <- map_dfr(1:100, read_s )
