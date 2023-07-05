library(tidyverse)
library(patchwork)
library(ggstance)
library(glue)
library(prismatic)
library(ggridges)
library(ggrastr)
library(here)

source(here("code/R/project_defaults.R"))

roh_length_threshold <- 1e3

import_genomes <- function(spec){
  read_tsv(here("data", "genomes", "filtered", glue("{spec}_filt.fa.gz.fai")),
           col_names = c("chr", "length", "offset", "linebases", "linewidth")) |> 
    mutate(spec = spec,
           idx = row_number(),
           end_pos = cumsum(length),
           start_pos = lag(end_pos, default = 0),
           eo = idx %% 2)
}

genomes <- specs |> map_dfr(import_genomes)

g_starts <- dplyr::select(genomes, chr, ref = spec, start_pos) 

# import_roh <- function(spec, partition){
#   read_tsv(here("results", "roh", "bcftools", glue("{spec}_{partition}_roh.tsv.gz")),
#            skip = 4,
#            col_names = c("RG", "sample", "chr", "start", "end",
#                          "Length_bp", "n_markers", "avg_phred_quality"))# |> 
#     # mutate(spec = spec,
#     #        sample = as.character(sample)) |> 
#     # left_join(g_starts, by = c("chr", "spec")) |> 
#     # mutate(g_start = start + start_pos,
#     #        g_end = end + start_pos)
# }

import_roh <- function(ref){
  read_tsv(here("results", "roh", "bcftools", "snp_based", glue("{ref}_roh.tsv.gz")),
           skip = 4,
           col_names = c("RG", "sample", "chr", "start", "end",
                         "Length_bp", "n_markers", "avg_phred_quality")) |> 
    mutate(ref = ref,
           sample = as.character(sample)) |> 
    left_join(g_starts, by = c("chr", "ref")) |> 
    mutate(g_start = start + start_pos,
           g_end = end + start_pos)
}
 |> |> |> |> |> |> |> |> 
data |> 
  filter(Length_bp <= 5)