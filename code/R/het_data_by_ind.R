library(tidyverse)
library(plyranges)
library(glue)
library(here)
library(data.table)
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
window_mirang <- slide_ranges(genome, width = window_width, step = window_step)

convert_012_to_hom <- \(mat){abs(1 - mat)}

run_window_het <- \(spec = "mirang", smpl_idx = 01){
  snp_ids <- read_tsv(glue("results/genotyping/012/mirang_filtered_{spec}_012.012.pos"),
                      col_names = c("chr", "pos")) |>
    mutate(snp = str_c(chr,"_",pos))
  inds <- read_lines(glue("results/genotyping/012/mirang_filtered_{spec}_012.012.indv"))
  
  geno_mat <- (fread(glue("results/genotyping/012/mirang_filtered_{spec}_012.tsv.gz"), sep = "\t")[smpl_idx,]) |> as.matrix()
  geno_mat[geno_mat == -1] <- NA
  
  granges_012 <- snp_ids |> 
    mutate(hom = convert_012_to_hom(geno_mat)[1,2:(dim(geno_mat)[2])],
           seqnames = str_remove(snp, "_[0-9]*$"),
           start = str_remove_all(snp, ".*_") |> as.integer()) |> 
    mutate(end = start) |> 
    select(-c(chr, snp)) |> 
    as_granges()
  
  avg_hom <- join_overlap_inner(window_mirang, granges_012) |> 
    as_tibble() |> 
    group_by(seqnames, start, end, strand, partition) |> 
    summarise(avg_hom = mean(hom, na.rm = TRUE),
              n_snps = n()) |> 
    mutate(mid = (start + end + 2) /2) |> 
    left_join(genome_tib |> dplyr::select(seqnames = chr, gstart = start_pos, eo)) |> 
    mutate(gpos = mid + gstart,
           ind = inds[smpl_idx],
           spec = spec)
  
  avg_hom
}

all_windows <- tibble(spec = rep(c("mirang", "mirleo"), c(20, 10)),
                      smpl_idx = c(1:20, 1:10)) |> 
  pmap_dfr(run_window_het)

write_tsv(all_windows,
          here(glue("results/het/win_het_ind_all_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.tsv.gz")))