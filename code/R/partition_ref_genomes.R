library(tidyverse)
library(patchwork)
library(prismatic)
library(ggforce)
library(geomtextpath)
library(glue)
library(here)

partition_genome <- \(spec, n_part = 20){
  # read in the scaffold sizes from genome index file
  data <- read_tsv(here("data", "genomes", glue("{spec}.fa.gz.fai")),
                   col_names = c("scaff", "length", "offset",
                                 "linebases", "linewidth"))
  
  # compute total genome length
  total_length <- sum(data$length)
  
  # compute genome wide end-coordinate for each scaffold,
  # assign to partition by integer division by 0.1 of
  # the total length and export as reference tsv file.  
  data2 <- data |>
    mutate(end = cumsum(length),
           partition_prep = str_pad(1 + end %/% (total_length/n_part + 1),
                                    width = 2,
                                    pad = 0),
           partition_idx = as.numeric(factor(partition_prep)))
  
  n_missing <- n_part - max(data2$partition_idx)
  # print(glue("n_missing: {n_missing}, n_part: {n_part}, max_part: {max(data2$partition_idx)}"))
  if(n_missing > 0){
    data3 <- data2 |> 
      group_by(partition_idx) |> 
      mutate(partition_adjust = if_else(partition_idx == max(partition_idx),
                                        partition_idx + (row_number() - 1) %/% ( ceiling(n() / (n_missing + 1) )),
                                        partition_idx)) |> 
      ungroup()
  } else {
    data3 <- data2 |> 
      mutate(partition_adjust = partition_idx)
  }
  
  data_out <- data3 |> 
    mutate(partition = str_pad(partition_adjust, width = 2, pad = 0)) |>
    dplyr::select(scaff, partition) |>
    write_tsv(file = here("data", "genomes", glue("{spec}_partitions.tsv")),
              col_names = FALSE)
}

check_paritions <- \(spec = "mirang"){
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
  
  p1 <- data |> 
    ggplot() +
    geom_vline(data = partitions,
               aes(xintercept = end_pos* 1e-9, color = partition),
               alpha =.3) +
    geom_step(aes(group = 1L, x = end_pos * 1e-9, y = idx, color = partition),
              linewidth = .3, alpha = .3) +
    geom_linerange(aes(group = 1L, xmin = lag(end_pos * 1e-9, default = 0),
                       xmax = end_pos * 1e-9,
                       y = idx, color = partition)) +
    labs(subtitle = spec)
  
  p2 <- partitions |> 
    ggplot(aes(x = partition, y = length, color = partition)) +
    geom_bar(stat = "identity", aes(fill = after_scale(clr_alpha(color)))) +
    geom_text(aes(y = length /2, label = n), angle = 90, hjust = -.2,  color = 'black')+
    labs(subtitle = "sequence within partitions")
  
  p3 <- partitions |> 
    ggplot(aes(x = partition, y = n, color = partition)) +
    geom_bar(stat = "identity", aes(fill = after_scale(clr_alpha(color)))) 

  p4 <- data |> 
    ggplot(aes(x = log10(length))) +
    geomtextpath::geom_textvline(data = tibble(x = c(100, 250, 500, 1000)),
                                 aes(xintercept = log10(x), 
                                     label = glue("{x} bp")),
               linetype = 3, color ="#006D2C", 
               size = 2.5, hjust = .85) +
    geom_density(color = "#C7E9C0",
                 fill = clr_alpha("#C7E9C0")) +
    scale_x_reverse() +
    coord_flip(expand = 0, xlim = log10(c(1e8, 30)))+
    labs(subtitle = "scaffold sizes")
  
  clr_ln <- 'gray50'
  p1 + (p2 / p3) + p4 +
    plot_layout(widths = c(.7, 1, .3)) +
    plot_annotation(title = spec) &
    theme_minimal() &
    theme(legend.position = "none",
          plot.subtitle = element_text(hjust = .5),
          panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(size = .3, color = clr_ln),
          axis.ticks.length = unit(3, 'pt'),
          axis.line = element_line(color = clr_ln, size = .4))
  }

n_part <- 20
c("mirang", "mirleo", 'filtered/mirang_filt', 'filtered/mirleo_filt') |>
  walk(partition_genome, n_part = n_part)

pp <- check_paritions('mirang') /
  check_paritions('mirleo') /
  check_paritions('filtered/mirang_filt') /
  check_paritions('filtered/mirleo_filt')

p <- pp &
 scale_color_manual(values = c(
   scales::colour_ramp(RColorBrewer::brewer.pal(6,"Greens")[6:2])(seq(0,1,length.out = ceiling(n_part/2))),
   scales::colour_ramp(RColorBrewer::brewer.pal(6,"Reds")[2:6])(seq(0,1,length.out = ceiling(n_part/2)))))

ggsave(filename = here("results", "img", "qc", "genome_partitions.pdf"),
       plot = p,
       width = 16,
       height = 16,
       device = cairo_pdf)