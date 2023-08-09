library(tidyverse)
library(prismatic)
library(patchwork)
library(glue)
library(here)
source(here("code/R/project_defaults.R"))

import_genomes <- function(spec){
  read_tsv(here("data", "genomes", "filtered", glue("{spec}_filt.fa.gz.fai")),
           col_names = c("chr", "length", "offset", "linebases", "linewidth")) |> 
    mutate(spec = spec,
           idx = row_number(),
           end_pos = cumsum(length),
           start_pos = lag(end_pos, default = 0),
           eo = idx %% 2)
}
import_pop <- \(spec){
  tibble(sample = read_lines(here(glue("results/pop/inds_{spec}.pop"))),
         spec = spec)}

samples <- c("mirang", "mirleo") |> 
  map_dfr(import_pop)


import_genomes <- function(spec){
  read_tsv(here("data", "genomes", "filtered", glue("{spec}_filt.fa.gz.fai")),
           col_names = c("chr", "length", "offset", "linebases", "linewidth")) |> 
    mutate(spec = spec,
           idx = row_number(),
           end_pos = cumsum(length),
           start_pos = lag(end_pos, default = 0),
           eo = idx %% 2)
}

genome <- import_genomes("mirang")

data_pheno <- read_tsv(here("results/pop/group_pheno_labeled.pop"),
                       col_names = c("sample_id", "phenotype"),
                       col_types = "cc") |> 
  mutate(phenotype = factor(phenotype, levels = c("worms", "control", "mirleo")))

data_het_ind <- read_tsv(here(glue("results/het/win_het_ind_all_w1Mb_s250kb.tsv.gz"))) |> 
  left_join(data_pheno, by = c(ind = "sample_id")) |> 
  rename(pheno = "phenotype") 


plink_path <- here("results/roh/plink/")
plink_folders <- dir(plink_path, pattern = "d[0-9]*$")


on_x <- read_tsv(here("results/genomes/sex_chrom/mirang_sex_chrom.bed"))
g_starts <- dplyr::select(genome, chr, ref = spec, start_pos, eo) 
# test <- "mirang_filtered_all_h0_wh3_n100_wn50_wm5_l10_g50_d50"
summarize_plink_run <- \(folder){
  
  pl_h <- str_replace(folder, ".*_h([0-9]*)_.*", "\\1")
  pl_wh <- str_replace(folder, ".*_wh([0-9]*)_.*", "\\1")
  pl_n <- str_replace(folder, ".*_n([0-9]*)_.*", "\\1")
  pl_wn <- str_replace(folder, ".*_wn([0-9]*)_.*", "\\1")
  pl_mis <- str_replace(folder, ".*_wm([0-9]*)_.*", "\\1")
  pl_l <- str_replace(folder, ".*_l([0-9]*)_.*", "\\1")
  pl_g <- str_replace(folder, ".*_g([0-9]*)_.*", "\\1")
  pl_den <- str_replace(folder, ".*_d([0-9]*)$", "\\1")
  
  plink_tag <- glue("h{pl_h}_wh{pl_wh}_n{pl_n}_wn{pl_wn}_wm{pl_mis}_l{pl_l}_g{pl_g}_d{pl_den}")
  
  plink_data <- read_table(glue("{plink_path}/mirang_filtered_all_{plink_tag}/mirang_filtered_all_{plink_tag}.hom"),
                           col_types = "ccdicciididdd") |> 
    mutate(CHR = str_remove(SNP1, ":.*"),
           autosome = !(CHR %in% on_x$chr)) |> 
    left_join(g_starts |> filter(ref == "mirang"), by = c(CHR = "chr")) |> 
    left_join(samples, by = c(IID = "sample")) |> 
    mutate(gstart = POS1 + start_pos,
           gend = POS2 + start_pos,
           gmid = (gstart + gend) / 2 )
  
  sum_bp <- genome |> 
    pluck("length") |> 
    sum()
  
  roh_summary_by_sample <- plink_data |>
    group_by(IID, spec) |> 
    summarize(sum_roh_length = sum(KB * 1e3)) |> 
    # left_join(mappable_lengths |> filter(roh_type == roh_version)) |> 
    ungroup() |> 
    mutate(f_roh = sum_roh_length / sum_bp) # |>
  # filter(Length_bp > roh_length_threshold) |> 
  
  tibble(data_plink = list(plink_data), roh_summary = list(roh_summary_by_sample))
  
}

plot_plink_comparison <- \(folder){
  data_plink <- summarize_plink_run(folder)
  plot_ind <- \(focal_ind = "ES2551"){
    data_bcftools <- read_tsv(here(glue("results/roh/bcftools/bed/max_callable/roh_max_{focal_ind}_on_mirang.bed")), col_names = c("chr", "start", "end")) |> 
      left_join(g_starts) |> 
      mutate(gstart = start_pos + start, 
             gend = start_pos + end,
             gmid = (gstart + gend) / 2 )
    
    data_het_ind |> 
      filter(#seqnames == genome$chr[7],
        ind == focal_ind) |> 
      ggplot() +
      geom_rect(data = genome |> select(-spec),
                aes(xmin = start_pos,
                    xmax = end_pos,
                    ymin = -Inf,
                    ymax = Inf,
                    fill = as.character(eo)),
                color = "transparent") +
      geom_hline(yintercept = 0:1, color = rgb(0,0,0,.2), linewidth = .2) +
      geom_line(aes(x = gpos, y = 1 - avg_hom, group = str_c(seqnames, pheno, ind), color = pheno), linewidth = .25, alpha = .6) +
      # geom_point(aes(x = gpos, y = 1 - avg_hom, group = seqnames, color = pheno), size = .15, shape = 19, alpha = .2) +
      # geom_point(data = data_plink$data_plink[[1]] |> filter(IID == focal_ind),
      #            aes(x = gmid, y = -.3, fill = after_scale(clr_lighten(color))), size = .75, shape = 21, alpha= .3) +
      geom_linerange(data = data_plink$data_plink[[1]] |> filter(IID == focal_ind),
                     aes(xmin = gstart, xmax = gend, y = -.3), lwd = 1.5, alpha= .7) +
      # geom_point(data = data_bcftools,
      #            aes(x = gmid, y = -.1, fill = after_scale(clr_lighten(color))), size = .75, shape = 21, alpha= .3) +
      geom_linerange(data = data_bcftools,
                     aes(xmin = gstart, xmax = gend, y = -.1), lwd = 1.5, alpha= .7) +
      scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(breaks = c(-.3, -.1, 0, 1), labels = c("plink", "bcftools", 0, 1)) +
      scale_color_manual(values = clr_pheno, guide = "none") +
      coord_cartesian(ylim = c(-.4, 1)) +
      labs(y = focal_ind) +
      theme_minimal(base_family = fnt_sel, base_size = 9) +
      theme(axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            # axis.text.y = element_blank(),
            panel.grid = element_blank())
  }
  
  ps <- data_pheno$sample_id |> 
    map(plot_ind)
  p <- wrap_plots(ps, ncol = 1)
  ggsave(plot = p, filename = glue("results/img/roh/comp_het/{folder}.png"), width = 10, height = 22)
}

plink_folders <- dir(plink_path, pattern = "d[0-9]*$")

plink_folders[4:48] |>  walk(plot_plink_comparison)
