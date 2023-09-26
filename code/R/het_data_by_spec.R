library(tidyverse)
library(plyranges)
library(glue)
library(here)
library(data.table)
source(here("code/R/project_defaults_shared.R"))

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

window_mirang <- slide_ranges(genome, width = window_width, step = window_step)

convert_012_to_hom <- \(mat){abs(1 - mat)}
convert_012_to_het <- \(mat){1 - abs(1 - mat)}

avg_het_by_snp <- \(spec){

  geno_mat <- (fread(glue("results/genotyping/012/mirang_filtered_{spec}_012.tsv.gz"), sep = "\t"))[,-1] |> as.matrix()
  geno_mat[geno_mat == -1] <- NA
  cum_het <- colSums(convert_012_to_het(geno_mat), na.rm = TRUE)
  n_typed <- colSums(!is.na(geno_mat), na.rm = TRUE)
  avg_het <- cum_het / n_typed
  avg_het
}

avg_het <- list(
  mirang = avg_het_by_snp("mirang"),
  mirleo = avg_het_by_snp("mirleo")
)

win_avg_het <- \(spec){
  snp_ids <- read_tsv(glue("results/genotyping/012/mirang_filtered_{spec}_012.012.pos"),
                      col_names = c("chr", "pos")) |>
    mutate(snp = str_c(chr,"_",pos))
  
  granges_012 <- snp_ids |> 
    mutate(mean_het = avg_het[[spec]],
           seqnames = str_remove(snp, "_[0-9]*$"),
           start = str_remove_all(snp, ".*_") |> as.integer()) |> 
    mutate(end = start) |> 
    select(-c(chr, snp)) |> 
    as_granges()
  
  join_overlap_inner(window_mirang, granges_012) |>
    as_tibble() |>
    group_by(seqnames, start, end, strand, partition) |>
    summarise(avg_het = mean(mean_het, na.rm = TRUE),
              win_het = sum(mean_het, na.rm = TRUE) / (end[1] - start[1])) |>
    mutate(mid = (start + end + 2) /2) |>
    left_join(genome_tib |>
                dplyr::select(seqnames = chr, gstart = start_pos, eo)) |>
    mutate(gpos = mid + gstart,
           spec = spec) |>
    ungroup()
}

windowed_avg_het <- c("mirang", "mirleo") |> 
  map_dfr(win_avg_het)

windowed_avg_het |> 
  write_tsv(here(glue("results/het/win_het_by_spec_w{window_width*1e-6}Mb_s{window_step*1e-3}kb.tsv.gz")))