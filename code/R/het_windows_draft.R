library(tidyverse)
library(plyranges)
library(glue)
library(here)
library(prismatic)
library(data.table)
library(ggrastr)
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
  
  p <- avg_hom |>
    ggplot() +
    geom_rect(data = genome_tib,
              aes(xmin = start_pos,
                  xmax = end_pos,
                  ymin = -Inf,
                  ymax = Inf,
                  fill = as.character(eo)),
              color = "transparent") +
    geom_line(aes(x = gpos, y = avg_hom, group = seqnames, color = spec), linewidth = .25) +
    geom_point(aes(x = gpos, y = avg_hom, group = seqnames, color = spec), size = .15, shape = 19, alpha = .2) +
    scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
    scale_x_continuous(expand = c(0, 0)) +
    # coord_cartesian(xlim = c(0, max(genome_tib$end_pos)),
    #                 ylim = c(-.1 , 1.1), 
    #                 expand = 0) +
    scale_color_manual(values = clrs, guide = "none") +
    labs(subtitle = glue("{inds[smpl_idx]} ({spec})")) +
    theme_minimal(base_family = fnt_sel, base_size = 9) +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          # axis.text.y = element_blank(),
          panel.grid = element_blank())

  ind_lab <- str_pad(smpl_idx, width = 2, pad = 0)
  
  ggsave(filename = here(glue("results/img/het/inds/win_n_snps_{spec}_{inds[smpl_idx]}_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.png")),
         plot = p, 
         width = 10, height = 2, bg = "white")
  
  avg_hom
}

# all_windows <- read_tsv(here(glue("results/het/win_het_ind_all_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.tsv.gz")))

all_windows <- tibble(spec = rep(c("mirang", "mirleo"), c(20, 10)),
       smpl_idx = c(1:20, 1:10)) |> 
  pmap_dfr(run_window_het)

all_windows |> 
  group_by(spec) |> 
  filter(as.numeric(as.factor(ind)) == 1 ) |> 
  ungroup() |> 
  ggplot() +
  geom_rect(data = genome_tib |> select(-spec), 
            aes(xmin = start_pos,
                xmax = end_pos, 
                ymin = -Inf,
                ymax = Inf, 
                fill = as.character(eo)),
            color = "transparent") +
  geom_line(aes(x = gpos, y = n_snps, group = seqnames, color = spec), linewidth = .25) +
  geom_point(aes(x = gpos, y = n_snps, group = seqnames, color = spec), size = .15, shape = 19, alpha = .2) +
  # scale_y_continuous(breaks = 0:1) +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') + 
  coord_cartesian(xlim = c(0, max(genome_tib$end_pos)),
                  # ylim = c(-.1 , 1.1),
                  expand = 0) +
  scale_color_manual(values = clrs, guide = "none") + 
  scale_x_continuous(expand = c(0, 0), labels = \(br){sprintf("%.1fGb", br * 1e-9)},
                     sec.axis = sec_axis(trans = identity, breaks = genome_tib$mid_pos[1:17], labels = 1:17)) +
  labs(subtitle = glue("snp density (w{sprintf('%.0fMb', window_width /1e6)} s{sprintf('%.0fkb', window_step /1e3)})")) +
  facet_grid(spec ~ ., switch = "y", scales = "free") +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.ticks.x.bottom = element_line(linewidth = .2),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 0))

ggsave(filename = here(glue("results/img/het/win_n_snps_ind_all_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.png")),
       # plot = p,
       width = 10, height = 4, bg = "white")

all_windows |> 
  ggplot() +
  geom_rect(data = genome_tib, 
            aes(xmin = start_pos,
                xmax = end_pos, 
                ymin = -Inf,
                ymax = Inf, 
                fill = as.character(eo)),
            color = "transparent") +
  geom_hline(yintercept = 0:1, color = rgb(0,0,0,.1), linewidth = .2) +
  geom_line(aes(x = gpos, y = avg_hom, group = seqnames, color = spec), linewidth = .25)  +
  geom_point(aes(x = gpos, y = avg_hom, group = seqnames, color = spec), size = .15, shape = 19, alpha = .2) +
  scale_y_continuous(breaks = 0:1) +
  scale_x_continuous(expand = c(0, 0), labels = \(br){sprintf("%.1fGb", br * 1e-9)},
                     sec.axis = sec_axis(trans = identity, breaks = genome_tib$mid_pos[1:17], labels = 1:17)) +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') + 
  coord_cartesian(xlim = c(0, max(genome_tib$end_pos)),
                  ylim = c(-.1 , 1.1), expand = 0) +
  scale_color_manual(values = clrs, guide = "none") + 
  # labs(subtitle = glue("{inds[smpl_idx]} ({spec})")) +
  facet_grid(ind ~ ., switch = "y") +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.ticks.x.bottom = element_line(linewidth = .2),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 0))

ggsave(filename = here(glue("results/img/het/win_hom_ind_all_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.png")),
       # plot = p,
       width = 10, height = 15, bg = "white")

write_tsv(all_windows,
          here(glue("results/het/win_het_ind_all_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.tsv.gz")))