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

# import_roh <- function(spec, partition){
#   read_tsv(here("results", "roh", "bcftools", glue("{spec}_{partition}_roh.tsv.gz")),
#            skip = 4,
#            col_names = c("RG", "sample", "chr", "start", "end",
#                          "Length_bp", "n_markers", "avg_phred_quality"))# |> 
#     # mutate(spec = spec,
#     #        sample = as.character(sample)) |> 
#     # left_join(g_starts, by = c("chr", "spec")) |> 
#     # mutate(g_start = start + start_pos,
#     #        g_end = end + start_pos)
# }

import_roh <- function(ref){
  read_tsv(here("results", "roh", "bcftools", "snp_based", glue("{ref}_roh.tsv.gz")),
           skip = 4,
           col_names = c("RG", "sample", "chr", "start", "end",
                         "Length_bp", "n_markers", "avg_phred_quality")) |> 
    mutate(ref = ref,
           sample = as.character(sample)) |> 
    left_join(g_starts, by = c("chr", "ref")) |> 
    mutate(g_start = start + start_pos,
           g_end = end + start_pos)
}

import_pop <- \(spec){
  tibble(sample = read_lines(here(glue("results/pop/inds_{spec}.pop"))),
         spec = spec)}

on_x <- read_tsv(here("results/genomes/sex_chrom/mirang_sex_chrom.bed"))
data <- import_roh("mirang") |> 
  left_join(map_dfr(c("mirang", "mirleo"), import_pop)) |> 
  mutate(on_x = chr %in% on_x$chr)
# import_roh("mirang","10")
# import_roh("mirang","20")
# 
# data <- str_pad(1:20, width = 2, pad = 0) |> 
#   map_dfr(import_roh, spec = "mirang")

#  >>>>>>>>>>> waiting for re-genotyping <<<<<<<<<<<<<
# data |> 
#   filter(chr %in% (c("NW_025578508.1", "NW_025578572.1", "NW_025578573.1") |> str_replace("NW_","JAAMPH"))) 

# ns <- read_tsv("ns.bed", col_names = c("chr", "start", "end")) |> 
#   mutate(length = end - start)
# 
# ns |> 
#   ggplot(aes(x = length)) +
#   geom_histogram(color = "black", fill = rgb(0,0,0,.2), bins = 35, boundary = 1)
genome_mirang <- genomes |> filter(spec == "mirang") |> mutate(mid_pos = (start_pos + end_pos) / 2)

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
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
  theme(panel.grid = element_blank())

data_summary <- data |> 
  filter(!on_x) |> 
  filter(Length_bp > roh_length_threshold) |> 
  group_by(sample, spec) |> 
  summarise(total_roh_length_mb = sum(Length_bp) * 1e-6,
            total_n_snp_in_roh = sum(n_markers) * 1e-3)# |> 
# left_join(genome)
# ,
#             f_roh = total_roh_length/20e9)

p2 <- data_summary |> 
  pivot_longer(cols = -c(sample, spec)) |> 
  ggplot(aes(y = 0, x = value, color = spec)) +
  geom_barh(stat = "identity",
            aes(fill = after_scale(clr_alpha(color)))) +
  facet_grid(sample ~ name, scale = "free", switch = "y")

p3 <- data |>
  filter(Length_bp > roh_length_threshold) |>
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
  scale_color_manual(values = clrs)
# scale_fill_manual(values = clrs |> clr_alpha() |> set_names(nm = names(clrs)))

p_out <- p1  + p3 + p2 + 
  plot_layout(widths = c(1, .6, .6), guides = "collect") +
  plot_annotation(title = "ROH sumary (max callable)",
                  caption = glue("only ROH larger {sprintf('%.0f',roh_length_threshold)} bp are included")) & 
  scale_color_manual(values = clrs) &
  theme_minimal(base_size = 9) &
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        strip.text.y.left = element_text(angle = 0),
        plot.caption = element_text(),
        plot.subtitle = element_text())

ggsave(filename = here("results/img/roh/mirang_roh_snps.pdf"),
       plot = p_out,
       width = 16, height = 6, device = cairo_pdf)

pq1 <- data |> 
  filter(Length_bp > roh_length_threshold) |> 
  ggplot(aes(x = avg_phred_quality)) +
  geom_vline(xintercept = 30, linetype = 3, color = rgb(0,0,0, .7)) +
  geom_density(aes(color = spec, fill = after_scale(clr_alpha(color)))) + 
  facet_wrap(sample ~ .) +
  scale_color_manual(values = clrs) +
  theme_minimal(base_size = 9)

pq2 <- data |> 
  filter(Length_bp > roh_length_threshold) |> 
  mutate(length_mb = Length_bp * 1e-6) |> 
  ggplot(aes(x = length_mb, y = avg_phred_quality)) +
  geom_hline(yintercept = 30, linetype = 3, color = rgb(0,0,0, .7)) +
  ggrastr::rasterise(geom_point(aes(color = spec), size = .3, alpha = .4)) + 
  facet_wrap(sample ~ .) +
  scale_color_manual(values = clrs) +
  theme_minimal(base_size = 9)

pq <- pq1 + pq2 + 
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = here("results/img/roh/mirang_roh_snps_qc.pdf"),
       plot = pq,
       width = 12, height = 5, device = cairo_pdf)
