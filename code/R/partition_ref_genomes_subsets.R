#!/usr/bin/env Rscript
# Rscript <spec> <n_part> <n_sub> <out_plt>
# Rscript filtered/mirang_filt 20 10 results/img/qc/partition_sub_{spec}.pdf
library(tidyverse)
library(patchwork)
library(glue)
library(here)

args <- commandArgs(trailingOnly = TRUE)
spec_base <- as.character(args[1])
spec <- str_c("filtered/", spec_base, "_filt")
n_part <- as.integer(args[2])
n_sub <- as.integer(args[3])
out_plt <- as.character(args[4])

partitions_base <- \(spec = "mirang"){
  data <- read_tsv(here("data", "genomes", glue("{spec}.fa.gz.fai")),
                   col_names = c("seq_name", "length", "offset", "linebases", "linewidth")) |> 
    left_join(read_tsv(here("data", "genomes", glue("{spec}_partitions.tsv")),
                       col_names = c("seq_name", 'partition'))) |> 
    mutate(idx = row_number(),
           end_pos = cumsum(length))
  
  partitions <- data |> 
    group_by(partition) |> 
    summarize(
      n = n(),
      length = sum(length),
      end_pos = max(end_pos),
      idx_start = min(idx),
      idx_end = max(idx)
    ) |> 
    ungroup() |> 
    mutate(failed = partition %in% c('02', '04', '06'))
  
  tibble(spec = spec, data = list(data), partitions = list(partitions))
}

sub_partitions <- \(data, part, n_subs = 10){
    partition_scafs <- data |> 
      filter(partition == part)
    
    n_scaffs <- nrow(partition_scafs)
    partition_cum_length <- sum(partition_scafs$length)
    
    sub_prep <- partition_scafs |> 
      mutate(share_subs = round(length / partition_cum_length * n_subs))
    
    subs_left <- n_subs - sum(sub_prep$share_subs)
    length_in_short_scaffolds <- sum(sub_prep$length[sub_prep$share_subs == 0])
    if(subs_left < 2){ 
      subs_break <- c(-Inf, Inf)
    } else {
      subs_break <- seq(0, length_in_short_scaffolds, length.out = subs_left + 1) 
    }
    
    sub_bins <- sub_prep |>
      mutate(is_0 = share_subs == 0) |>
      group_by(is_0) |>
      mutate(bin = as.integer((cut(cumsum(length), breaks = subs_break)))) |> 
      ungroup() |> 
      mutate(cum_subs = cumsum(share_subs),
             cum_subs_bin = if_else(is_0, cum_subs + bin, cum_subs)) |> 
      group_by(cum_subs_bin) |> 
      mutate(n_scaffolds = length(seq_name)) |> 
      ungroup()
    
    # make sure that the total number of partition subs matches the specified n
    sub_sum_check <- max(sub_bins$cum_subs_bin) 
    # sub_sum_check
    if( sub_sum_check != n_subs){stop(glue("n sub partitions does not match for partition {part}! ({sub_sum_check} instead of {n_subs})"))}

    make_breaks <- \(n_scaffolds, length, share_subs,...){
      if(n_scaffolds > 1){
        tibble(start = 1, end = length)
      } else {
         tibble(start = ceiling(seq(0, length, length.out = share_subs+1))[1:share_subs] + 1,
                 end = ceiling(seq(0, length, length.out = share_subs+1))[2:(share_subs+1)])
        }
    }

    sub_bins |>
      mutate(breaks = pmap(sub_bins, make_breaks)) |>
      unnest(breaks) |>
      mutate(sub_idx = if_else(is_0, cum_subs_bin, row_number()),
             sub_length = end - start + 1) |>
      select(seq_name, sub_idx, start, end, sub_length) |>
      mutate(part = part) |>
      select(part, everything())
  }

export_intervalls <- \(data, part_q, sub, spec){
  dir.create(as.character(here("data", "genomes", glue("{spec}_partitions"))),
             showWarnings = FALSE)
  data |> 
    filter(part == part_q, sub_idx == sub) |> 
    mutate(interv = str_c(seq_name,":",start,"-",end)) |> 
    select(interv) |> 
    write_tsv(file = here("data", "genomes", glue("{spec}_partitions"), glue("part_{part_q}_sub_{sub}.intervals")),
              col_names = FALSE)
}

part_labs <- str_pad(1:n_part, width = 2, pad = "0")
data <- partitions_base(spec)$data[[1]]

all_subs <- part_labs |> 
    map_dfr(sub_partitions, n_subs = n_sub, data = data)

crossing(part_q = part_labs,
         sub = 1:n_sub) |> 
  pwalk(export_intervalls,
        data = all_subs,
        spec = spec)
 
subs_summary <- all_subs |>
  group_by(part, sub_idx) |> 
  summarise(length = sum(sub_length),
            n_scaff = length(seq_name))

p1 <- subs_summary |> 
  ggplot(aes(x = sub_idx, y = length)) +
  geom_bar(stat = 'identity', color = "gray60", fill = "gray90") +
  facet_wrap(part ~ ., scale = "free")

p2 <- subs_summary |> 
  ggplot(aes(x = sub_idx, y = n_scaff)) +
  geom_bar(stat = 'identity', color = "gray60", fill = "gray90") +
  facet_wrap(part ~ ., scale = "free")

p_out <- p1 + p2  +
  plot_annotation(title = glue("{spec} partition subs")) &
  theme_minimal()

ggsave(filename = out_plt,
       plot = p_out,
       width = 15,
       height = 7,
       device = cairo_pdf)