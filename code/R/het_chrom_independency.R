library(tidyverse)
library(plyranges)
library(glue)
library(here)
library(prismatic)
library(GGally)
source(here("code/R/project_defaults.R"))

import_genomes <- function(spec){
  read_tsv(here("data", "genomes", "filtered", glue("{spec}_filt.fa.gz.fai")),
           col_names = c("chr", "length", "offset", "linebases", "linewidth")) |> 
    mutate(spec = spec,
           idx = row_number(),
           end_pos = cumsum(length),
           start_pos = lag(end_pos, default = 0),
           mid_pos = (end_pos + start_pos)/2,
           eo = idx %% 2)
}

genome_tib <- import_genomes("mirang")
genome <- genome_tib  |> 
  mutate(start = 0) |> 
  dplyr::select(seqnames = chr, start, end = length) |> 
  as_granges()

window_width <- 1e6
window_step <- 25e4

all_windows <- read_tsv(here(glue("results/het/win_het_ind_all_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.tsv.gz"))) |> 
  filter(start %% 1e6 == 0) |> 
  mutate(length = end - start + 1,
         w_het = avg_hom * n_snps,
         win_het = w_het / length) |> 
  group_by(seqnames, ind, spec) |> 
  summarise(chr_het = weighted.mean(win_het, length))

all_windows |> 
  filter(spec == "mirang") |> 
  select(-spec) |>
  left_join(genome_tib |> select(seqnames = chr, length)) |> 
  arrange(-length) |> 
  ungroup() |> 
  filter(as.numeric(as.factor(-length)) < 18) |> 
  select(-length) |> 
  pivot_wider(names_from = seqnames,
              values_from = chr_het)  |> 
  select(-ind) |> 
  ggpairs()


ggsave("results/img/het/chr_independency.pdf", width = 16, height = 14, device = cairo_pdf)
