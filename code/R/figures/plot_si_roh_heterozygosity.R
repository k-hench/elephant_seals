library(tidyverse)
library(prismatic)
library(glue)
library(here)
library(ggrastr)
source(here("code/R/project_defaults.R"))

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
  rename(pheno = "phenotype") |> 
  mutate(ind = fct_reorder(ind, as.numeric(factor(str_c(as.numeric(pheno),"-",ind)))))

on_x <- read_tsv(here("results/genomes/sex_chrom/mirang_sex_chrom.bed"))
g_starts <- dplyr::select(genome, chr, ref = spec, start_pos, eo) 

inds <- read_tsv(here("data/file_info.tsv")) |>
  filter(!duplicated(sample_id)) |>
  pluck("sample_id")

read_bcftools <- \(focal_ind){
  read_tsv(here(glue("results/roh/bcftools/bed/max_callable/roh_max_{focal_ind}_on_mirang.bed")), col_names = c("chr", "start", "end")) |> 
    left_join(g_starts) |> 
    mutate(gstart = start_pos + start, 
           gend = start_pos + end,
           gmid = (gstart + gend) / 2,
           ind = focal_ind)
}

data_bcftools <- map_dfr(inds, read_bcftools) |> 
  mutate(ind = factor(ind, levels = levels(data_het_ind$ind)))
  
p <- data_het_ind |> 
      ggplot() +
      geom_vline(data = genome |> select(-spec) |> head(20),
                 aes(xintercept = end_pos),
                 alpha = .1) +
      geom_hline(yintercept = 0:1,
                 color = rgb(0,0,0,.2), linewidth = .2) +
      facet_wrap(ind ~., ncol = 2, strip.position = "right") +
      rasterise(geom_line(aes(x = gpos, y = 1 - avg_hom,
                    group = str_c(seqnames, pheno, ind), color = pheno),
                    linewidth = .2, alpha = .6), dpi = 300) +
      rasterise(geom_linerange(data = data_bcftools |> 
                                 filter((end - start) > 1e6),
                     aes(xmin = gstart, xmax = gend, y = -.2), lwd = 1.5, alpha= .7), dpi = 300) +
      scale_x_continuous("Position on genome", expand = c(0, 0), labels = \(x){sprintf("%.1f Gb",x *1e-9)}) +
      scale_y_continuous(breaks = c(0, 1), labels = c(0, 1)) +
      scale_color_manual("Phenotype", values = clr_pheno) +
      coord_cartesian(ylim = c(-.4, 1)) +
      labs(y = "Heterozygosity within variable sites") +
      theme_ms() +
      theme(strip.background = element_blank(),
            legend.position = "bottom",
            strip.text = element_text(size = 6))

ggsave(filename = here("results/img/final/sf_het.pdf"),
       plot = p,
       width = 12, height = 7,
       device = cairo_pdf)

ggsave(filename = here("results/img/final/sf_het.png"),
       plot = p,
       width = 12, height = 7)
