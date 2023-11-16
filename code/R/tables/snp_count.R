library(tidyverse)
library(here)
library(ggstance)
library(prismatic)
library(glue)
source("code/R/project_defaults.R")

ref <- "mirang"

files <- c(here(glue("results/genotyping/raw/{ref}_raw_snps_stat.txt")),
           dir(here("results/genotyping/filtered/"),
               pattern = glue("{ref}.*.txt"), full.names = TRUE))

get_vcfstats <- \(file){
  tibble(raw = read_lines(file) |> 
           str_remove("kept ") |> 
           str_replace(" out of [a-z ]*", ",") |> 
           str_replace(" ", ",")) |> 
    separate(raw, into = c("n","n_total", "stat"), convert = TRUE) |> 
    select(-n) |> 
    mutate(file = str_remove(file, ".*/") |> str_remove("_stat.txt")) |> 
    pivot_wider(names_from = stat,values_from = "n_total")
}

files |> 
  map_dfr(get_vcfstats) |> 
  mutate(file = factor(file, levels = c("mirang_raw_snps", "mirang_filtered_all",
                                        "mirang_bi-allelic", "mirang_mac1",
                                        "mirang_filtered_mirang", "mirang_filtered_mirleo") |> 
                         rev())) |> 
  arrange(-as.numeric(file))
