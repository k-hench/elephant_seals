library(tidyverse)
library(glue)
library(prismatic)
library(patchwork)
library(here)
source(here("code/R/project_defaults.R"))

import_genomes <- function(spec){
  read_tsv(here("data", "genomes", "filtered", glue("{spec}_filt.fa.gz.fai")),
           col_names = c("chr", "length", "offset", "linebases", "linewidth")) |> 
    mutate(spec = spec,
           idx = row_number(),
           end_pos = cumsum(length),
           start_pos = lag(end_pos, default = 0),
           mid_pos = (start_pos + end_pos) / 2,
           eo = idx %% 2)
}

genome <- import_genomes("mirang")

data_pheno <- read_tsv("results/pop/group_pheno_labeled.pop",
                       col_names = c("sample_id", "phenotype"),
                       col_types = "cc") |> 
  mutate(phenotype = factor(phenotype, levels = c("worms", "control", "mirleo")))

read_pi_spec <- \(part){
  part_pad <- str_pad(part, width = 2, pad = 0)
  read_csv(here(glue("results/pi/mirang_pi_dxy_{part_pad}.tsv.gz")))
}

read_pi_pheno <- \(part){
  part_pad <- str_pad(part, width = 2, pad = 0)
  read_csv(here(glue("results/pi/mirang_pheno_pi_dxy_{part_pad}.tsv.gz")))
}

data_to_genome <- \(data){
  data |> 
    left_join(genome |> select(scaffold = chr, start_pos, eo)) |> 
    mutate(gstart = start + start_pos,
           gend = end + start_pos,
           gmid = (gstart + gend)/2) |> 
    select(scaffold, gmid, sites, starts_with("pi")) |> 
    pivot_longer(cols = starts_with("pi"),
                 names_to = "pheno",
                 names_transform = \(str){str_remove(str, "pi_")},
                 values_to = "pi")
}

data_pi_spec <- map_dfr(1:20, read_pi_spec) |> data_to_genome()
data_pi_pheno <- map_dfr(1:20, read_pi_pheno) |> data_to_genome() |> 
  mutate(pheno = factor(pheno, levels = c("worms", "control", "mirleo")))

p1 <- ggplot() +
  geom_rect(data = genome |> select(-spec) |>  filter(row_number() < 18),
            aes(xmin = start_pos,
                xmax = end_pos,
                ymin = -Inf,
                ymax = Inf,
                fill = as.character(eo)),
            color = "transparent") +
  geom_point(data = data_pi_spec |> filter(scaffold %in% genome$chr[1:17]),
             aes(x = gmid, y = pi, color = pheno), size = .05, shape = 19, alpha = .2) +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
  scale_color_manual(values = clrs, guide = "none") +
  scale_x_continuous(expand = c(0, 0), labels = \(br){sprintf("%.1fGb", br * 1e-9)},
                     sec.axis = sec_axis(trans = identity, breaks = genome$mid_pos[1:17], labels = 1:17)) +
  scale_y_continuous("average pi (\U03C0)") +
  facet_grid(pheno ~ ., switch = "y") +
  labs(subtitle = glue("\U03C0 by species (w100kb s25kb)")) +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(strip.placement = "outside",
        axis.ticks.y = element_blank(),
        axis.ticks.x.bottom = element_line(linewidth = .2),
        axis.title.x = element_blank(),
        panel.grid = element_blank())

ggsave(filename = here(glue("results/img/pi/win_pi_spec_w100kb_s25kb.png")),
       plot = p1,
       width = 10, height = 3, bg = "white")


p2 <- ggplot() +
  geom_rect(data = genome |> select(-spec) |>  filter(row_number() < 18),
            aes(xmin = start_pos,
                xmax = end_pos,
                ymin = -Inf,
                ymax = Inf,
                fill = as.character(eo)),
            color = "transparent") +
  geom_point(data = data_pi_pheno |> filter(scaffold %in% genome$chr[1:17]),
             aes(x = gmid, y = pi, color = pheno), size = .05, shape = 19, alpha = .2) +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
  scale_color_manual(values = clr_pheno, guide = "none") +
  scale_x_continuous(expand = c(0, 0), labels = \(br){sprintf("%.1fGb", br * 1e-9)},
                     sec.axis = sec_axis(trans = identity, breaks = genome$mid_pos[1:17], labels = 1:17)) +
  scale_y_continuous("average pi (\U03C0)") +
  facet_grid(pheno ~ ., switch = "y") +
  labs(subtitle = glue("\U03C0 by phenotype (w100kb s25kb)")) +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(strip.placement = "outside",
        axis.ticks.y = element_blank(),
        axis.ticks.x.bottom = element_line(linewidth = .2),
        axis.title.x = element_blank(),
        panel.grid = element_blank())

ggsave(filename = here(glue("results/img/pi/win_pi_pheno_w100kb_s25kb.png")),
       plot = p2,
       width = 10, height = 5, bg = "white")

