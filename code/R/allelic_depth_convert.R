#!/usr/bin/env Rscript
# Rscript <tsv_file> <out_file> <n_samples>
# Rscript mirang_filtered_mirang_ad.tsv.gz mirang_filtered_mirang_het 20
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
tsv_file <- as.character(args[1])
out_file <- as.character(args[2])
n_samples <- as.integer(args[3])

read_tsv(tsv_file,
         col_types = str_c(c("c", "d", rep("c", n_samples + 1)), collapse = "")) |> 
  pivot_longer(ends_with(".AD"),
               names_to = "ind",
               names_transform = \(str){str_remove(str, ".AD$")}) |>
  mutate(ref = as.integer(str_remove(value, ",.*")),
         alt = as.integer(str_remove(value, ".*,")),
         major = if_else(ref > alt, ref, alt),
         minor = if_else(ref > alt, alt, ref)) |> 
  group_by(ind, minor, major) |> 
  count() |> 
  ungroup() |> 
  filter(minor > 0) |>
  write_tsv(file =  out_file)
