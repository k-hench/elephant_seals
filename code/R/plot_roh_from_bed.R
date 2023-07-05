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
prefix <- c(certain = "cert", callable = "max")
roh_version <- "certain" # "callable" #

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
"results/roh/bcftools/bed/max_callable/"

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

# import_roh("ES2816", "mirang")

on_x <- read_tsv(here("results/genomes/sex_chrom/mirang_sex_chrom.bed"))
data <- samples$sample |>
  map_dfr(import_roh, ref = "mirang") |> 
  mutate(on_x = chr %in% on_x$chr)

genome_mirang <- genomes |> filter(spec == "mirang") |> mutate(mid_pos = (start_pos + end_pos) / 2)
# ggplot() +
#   geom_rect(data =genome_mirang , 
#             aes(xmin = start_pos,
#                 xmax = end_pos, 
#                 ymin = -Inf,
#                 ymax = Inf, 
#                 fill = as.character(eo)),
#             color = "transparent") +
#   scale_x_continuous(breaks = genome_mirang$mid_pos[1:17],
#                      label = genome_mirang$chr[1:17] |> str_remove("N[CW]_0723")) +
#   scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
#   coord_cartesian(ylim = c(-1, 1)) +
#   theme_minimal()

p1 <- data |> 
  filter(!on_x) |> 
  filter(Length_bp > roh_length_threshold) |> 
  ggplot() +
  geom_rect(data = genomes |> filter(spec == "mirang"), 
            aes(xmin = start_pos,
                xmax = end_pos, 
                ymin = -Inf,
                ymax = Inf, 
                fill = as.character(eo)),
            color = "transparent") +
  rasterise(geom_linerange(aes(xmin = g_start, xmax = g_end, y = 0, color = spec),
                           linewidth = 4, alpha = .7), dpi = 450) +
  rasterise(geom_linerange(data = data |> 
                             filter(on_x, Length_bp > roh_length_threshold),
                           aes(xmin = g_start, xmax = g_end, y = 0), color = "gray85",
                           linewidth = 4, alpha = .9), dpi = 450) +
  scale_x_continuous(breaks = genome_mirang$mid_pos[1:17],
                     label = genome_mirang$chr[1:17] |> str_remove("N[CW]_0723")) +
  # geom_linerange(data = ns, 
  #                aes(xmin = start, xmax = end, y = 0), linewidth = 2, alpha = 1, color = "red") +
  facet_grid(sample ~ . , switch = "y") +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none')+ 
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank())

data_summary <- data |> 
  filter(!on_x) |> 
  filter(Length_bp > roh_length_threshold) |> 
  group_by(sample, spec) |> 
  summarise(total_roh_length_mb = sum(Length_bp) * 1e-6)# |> 
# left_join(genome)
# ,
#             f_roh = total_roh_length/20e9)

p2 <- data_summary |> 
  pivot_longer(cols = -c(sample, spec)) |> 
  ggplot(aes(y = 0, x = value, color = spec)) +
  geom_barh(stat = "identity",
            aes(fill = after_scale(clr_alpha(color)))) +
  facet_grid(sample ~ name, scale = "free", switch = "y")+ 
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

p3 <- data |>
  filter(Length_bp > roh_length_threshold) |>
  filter(!on_x) |> 
  ggplot(aes(x = log10(Length_bp),
             y = sample,
             color = spec)) +
  # geom_boxploth(aes(color = spec,
  #                   fill = spec,
  #                   group = sample))+
  geom_density_ridges2(aes(fill = after_scale(clr_alpha(color)),
                           group = sample)) +
  # geom_density_ridges(stat = "binline", bins = 45, scale = 0.95, draw_baseline = FALSE,
  #                     aes(fill = after_scale(clr_alpha(color)),
  #                         group = sample))+
  facet_grid(sample ~ ., scale = "free", switch = "y") +
  scale_color_manual(values = clrs) +
  coord_cartesian(xlim = c(2.8, 7.3), expand = 0) + 
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
# scale_fill_manual(values = clrs |> clr_alpha() |> set_names(nm = names(clrs)))



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

p4 <- roh_summary_by_sample |> 
  ggplot(aes(x = spec, color = spec)) +
  ggbeeswarm::geom_beeswarm(aes(y = f_roh), alpha = .7) +
  scale_color_manual(values = clrs) + 
  theme_minimal(base_family = fnt_sel, base_size = 9)

p_out <- p1  + p3 + p2 + p4 + 
  plot_layout(widths = c(1, .6, .3, .3), nrow = 1, guides = "collect") +
  plot_annotation(title = glue("ROH sumary (max {roh_version})"),
                  caption = glue("only ROH larger {sprintf('%.0f',roh_length_threshold)} bp on autosomes are included")) & 
  scale_color_manual(values = clrs) &
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        strip.text.y.left = element_text(angle = 0),
        plot.caption = element_text(),
        plot.subtitle = element_text())

ggsave(here(glue("results/img/roh/mirang_roh_max_{roh_version}_bed.pdf")),
       width = 16, height = 6, device = cairo_pdf)

ref <- "mirang"
roh_summary_by_sample |> 
  rename(roh_length_bp = "sum_roh_length",
         reference_length_bp = "sum_bp") |> 
  select(sample, spec, ref, roh_type, roh_length_bp, reference_length_bp, f_roh) |> 
  write_tsv(here("results", "roh", "bcftools", glue("summary_roh_{roh_version}_by_sample_on_{ref}.tsv")))


