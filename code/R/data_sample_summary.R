library(tidyverse)
library(here)

data <- read_tsv(here("data/file_info.tsv"))

data |> 
  filter(!duplicated(sample_id)) |> 
  select(spec, sample_id, treatment) |> 
  pluck("treatment") |> 
  table()

data |> 
  filter(!duplicated(sample_id)) |> 
  select(spec, sample_id, treatment) |> 
  write_tsv("~/Dropbox/Elephant seal paper/1 Data/wgs_samples.tsv")
