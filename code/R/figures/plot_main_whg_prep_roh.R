library(tidyverse)
library(patchwork)
library(ggstance)
library(glue)
library(prismatic)
library(ggridges)
library(here)
source(here("code/R/project_defaults.R"))
source(here("code/R/project_defaults_shared.R"))

roh_length_threshold <- 1e6 #1e3
prefix <- c(certain = "cert", callable = "max")
roh_version <- "callable" # "certain" # 

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

import_pop <- \(spec){
  tibble(sample = read_lines(here(glue("results/pop/inds_{spec}.pop"))),
         spec = spec)}

samples <- c("mirang", "mirleo") |> 
  map_dfr(import_pop)

import_roh <- function(sample, ref){
  read_tsv(here("results", "roh", "bcftools", "bed", glue("max_{roh_version}"),glue("roh_{prefix[roh_version]}_{sample}_on_{ref}.bed")),
           col_names = c("chr", "start", "end")) |> 
    mutate(Length_bp = end - start,
           ref = ref,
           sample = as.character(sample)) |> 
    left_join(g_starts, by = c("chr", "ref")) |> 
    mutate(g_start = start + start_pos,
           g_end = end + start_pos) |> 
    left_join(samples)
}

on_x <- read_tsv(here("results/genomes/sex_chrom/mirang_sex_chrom.bed"))
clr_lab <- c(mirang = clr_default[[1]], mirleo = clr_default[[2]])

data <- samples$sample |>
  map_dfr(import_roh, ref = "mirang") |> 
  mutate(on_x = chr %in% on_x$chr,
         sample_ord = as.numeric(factor(str_c(spec,"_" , str_pad(str_remove(sample, "[A-Z]*"), width = 7, pad = "0")))),
         sample_lab = str_c(glue("<span style='color:{clr_lab[spec]}'>"),sample,"</span>"))

# numeric summary for reduced F_ROH estimates
genome_mirang <- genomes |>
  filter(spec == "mirang") |>
  mutate(mid_pos = (start_pos + end_pos) / 2)

mappable_genome_length <- \(sample, ref, roh_type){
  if(roh_type == "certain"){
    data <- read_tsv(here("results/qc/coverage/masks/",glue("{sample}_on_{ref}_binary_covmask.bed.gz")),
                     col_names = c("chr", "start", "end")) |> 
      mutate(Length_bp = end - start) |> 
      filter(!(chr %in% on_x$chr)) |> 
      summarise(sum_bp = sum(Length_bp))
  } else if(roh_type == "callable") {
    data <- genomes |> 
      filter(spec == ref) |>
      filter(!(chr %in% on_x$chr)) |> 
      summarise(sum_bp = sum(length))
  } 
  data |> 
    mutate(ref = ref,
           sample = as.character(sample),
           roh_type = roh_type) |> 
    left_join(samples)
}

mappable_lengths <- crossing(roh_type = c("certain", "callable"),
                             sample = samples$sample) |> 
  pmap_dfr(mappable_genome_length, ref = "mirang")

roh_summary_by_sample <- data |>
  filter(Length_bp > roh_length_threshold) |> 
  group_by(sample, ref) |> 
  summarize(sum_roh_length = sum(Length_bp)) |> 
  left_join(mappable_lengths |> filter(roh_type == roh_version)) |> 
  ungroup() |> 
  mutate(f_roh = sum_roh_length / sum_bp)

data_phenotype <- read_tsv(here("data/file_info.tsv")) |> 
  filter(!duplicated(sample_id)) |> 
  select(sample = sample_id, treatment) |> 
  mutate(treatment = replace_na(treatment, "mirleo"))

p1 <- roh_summary_by_sample |> 
  left_join(data_phenotype) |> 
  ggplot(aes(x = spec_names[spec], y = f_roh)) +
  geom_jitter(aes(color = treatment,
               fill = after_scale(clr_alpha(color))),
              size = point_sz,
              height = 0, shape = 21, width = .25) +
  scale_color_manual(values = c(clr_pheno, mirleo = clr_default[[2]]),
                     labels = c(lab_pheno, mirang = "Southern elephant seal"),
                     guide = "none") +
  labs(y = "*F*<sub>ROH</sub>",
       x = NULL)+
  theme_ms()+
  theme(axis.title.x = ggtext::element_markdown(family = fnt_sel),
        axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.text.x = element_text(face = "italic"))

saveRDS(object = p1,
        here(glue("results/img/R/p_f_rho_{roh_version}_ind.Rds")))

roh_summary_by_sample |> 
  group_by(spec) |> 
  summarise(mean_froh = mean(f_roh),
            sd_froh = sd(f_roh)) |> 
  mutate(across(-spec, \(x){sprintf("%.2f", x)}))

genome_length_autosomes <- genomes |> 
  filter(spec == "mirang",
         !(chr %in% on_x$chr)) |> 
  pluck("length") |> 
  sum()
  
plot_cum_froh <- \(roh_length_threshold = 1e3, n_classes = 301){
  data |>
    filter(Length_bp > roh_length_threshold) |>
    filter(!on_x) |> 
    group_by(sample) |> 
    arrange(sample, Length_bp) |> 
    mutate(length_class = cut(Length_bp, breaks = (10  ^ seq(3,9, length.out = n_classes + 1)),
                              right = FALSE)) |>
    group_by(spec, sample, length_class) |> 
    summarise(sum_roh = sum(Length_bp)) |> 
    ungroup() |> 
    mutate(length_class = str_remove(length_class, "\\[") |> str_remove("\\)")) |> 
    separate(length_class, into = c("min_length", "max_length"), sep = ",", convert = TRUE) |> 
    mutate(length_class = (log10(min_length) + log10(max_length)) / 2) |> 
    group_by(spec, sample) |> 
    mutate(cum_roh_class = cumsum(sum_roh),
           cum_froh = cum_roh_class/ genome_length_autosomes ) |> 
    ggplot(aes(x = length_class,
               #x = Length_bp,
               y = cum_froh,
               color = spec,
               group = sample)) +
    geom_line(alpha = .45) +
    scale_color_manual(values = c(mirang = clr_default[[1]], 
                                  mirleo = clr_default[[2]]),
                       guide = "none") +
    scale_x_continuous(breaks = 3:7,
                       labels = c("1kb", "10kb", "100kb", "1Mb", "10Mb")) +
    coord_cartesian(xlim = c(2.8, 7.3),
                    ylim = c(0, .63),
                    expand = 0) +
    labs(x = "ROH Length",
         y = "Cumulative *F*<sub>ROH</sub>",
         subtitle = glue("min ROH length: {sprintf('%.0f',roh_length_threshold*1e-3)}kb"))+
    theme_ms() +
    theme(axis.title.y = ggtext::element_markdown(family = fnt_sel),
          axis.title.x = ggtext::element_markdown(family = fnt_sel),
          plot.subtitle = element_text(hjust = .5))
}

p_cum_roh <- 10^(3:6) |> map(plot_cum_froh)
p_cum_out <- p_cum_roh |> wrap_plots(nrow = 1)

saveRDS(object = p_cum_out,
        here(glue("results/img/R/p_cum_f_rho_{roh_version}_thresholds.Rds")))

saveRDS(object = p_cum_roh[[1]] + labs(subtitle = NULL),
        here(glue("results/img/R/p_cum_f_rho_{roh_version}.Rds")))

p_cum_pheno <- data |>
  filter(Length_bp > roh_length_threshold) |>
  filter(!on_x) |> 
  group_by(sample) |> 
  arrange(sample, Length_bp) |> 
  mutate(length_class = cut(Length_bp, breaks = (10  ^ seq(3,9, length.out = 301 + 1)),
                            right = FALSE)) |>
  group_by(spec, sample, length_class) |> 
  summarise(sum_roh = sum(Length_bp)) |> 
  ungroup() |> 
  mutate(length_class = str_remove(length_class, "\\[") |> str_remove("\\)")) |> 
  separate(length_class, into = c("min_length", "max_length"), sep = ",", convert = TRUE) |> 
  mutate(length_class = (log10(min_length) + log10(max_length)) / 2) |> 
  group_by(spec, sample) |> 
  mutate(cum_roh_class = cumsum(sum_roh),
         cum_froh = cum_roh_class/ genome_length_autosomes ) |> 
  left_join(data_phenotype) |> 
  ggplot(aes(x = length_class,
             #x = Length_bp,
             y = cum_froh,
             color = treatment,
             group = sample)) +
  geom_line(alpha = .45) +
  scale_color_manual(values = c(clr_pheno, mirleo = clr_default[[2]]),
                     labels = c(lab_pheno, mirang = "Southern elephant seal"),
                     guide = "none") +
  scale_x_continuous(breaks = 3:7,
                     labels = c("1kb", "10kb", "100kb", "1Mb", "10Mb")) +
  coord_cartesian(xlim = c(2.8, 7.3),
                  ylim = c(0, .63),
                  expand = 0) +
  labs(x = "ROH length",
       y = "Cumulative *F*<sub>ROH</sub>",
       subtitle = glue("min ROH length: {sprintf('%.0f',roh_length_threshold*1e-3)}kb"))+
  theme_ms() +
  theme(axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.title.x = ggtext::element_markdown(family = fnt_sel),
        plot.subtitle = element_text(hjust = .5))

saveRDS(object = p_cum_pheno + labs(subtitle = NULL),
        here(glue("results/img/R/p_cum_f_rho_{roh_version}_pheno.Rds")))

# numeric summary for reduced ROH length threshold
data |>
  filter(Length_bp > 1e3) |> 
  group_by(sample, ref) |> 
  summarize(sum_roh_length = sum(Length_bp)) |> 
  left_join(mappable_lengths |> filter(roh_type == roh_version)) |> 
  ungroup() |> 
  mutate(f_roh = sum_roh_length / sum_bp) |> 
  group_by(spec) |> 
  summarise(mean_froh = mean(f_roh),
            sd_froh = sd(f_roh)) |> 
  mutate(across(-spec, \(x){sprintf("%.2f", x)}))
