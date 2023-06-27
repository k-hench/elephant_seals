#!/usr/bin/env Rscript
# Rscript <tsv_file> <out_plt> <max_n_ref>
# Rscript mirang_filtered_mirang_het_bin.tsv.gz mirang_filtered_mirang_het_stats.pdf 30
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
tsv_file <- as.character(args[1])
out_plt <- as.character(args[2])
max_n_ref <- as.integer(args[3])

library(tidyverse)
library(prismatic)
library(patchwork)

data <- read_tsv(tsv_file, col_types = "ciii")

p1 <- data |> 
  group_by(ind) |> 
  summarise(n = sum(n)) |> 
  ggplot(aes(x = n *1e-3)) +
  geom_histogram(boundary = 0, color = "gray70", fill = "gray95")+
  labs(subtitle = "n SNPs per sample", x = "n SNPs (k)", y = "sample count") +
  theme_light() +
  theme(plot.subtitle = element_text(hjust = .5))

p2 <- data |> 
  ggplot(aes(x = major, y = minor, color = n)) +
  geom_tile(aes(fill = after_scale(clr_lighten(color)))) +
  geom_abline(data = tibble(a = c(0.5, 0.3, 0.2)),
              aes(slope = a*2, intercept = 0, linetype = factor(a)),
              alpha = .5, linewidth = .4) +
  facet_wrap(ind ~ .) +
  scale_color_viridis_c() +
  scale_linetype_manual(values = c(`0.5` = 1, `0.3` = 2, `0.2` = 3),
                        guide = "none")+
  coord_equal(xlim = c(0, max_n_ref),
              ylim = c(0, max_n_ref * .75)) +
  labs(subtitle = "distribution of allele frequency", x = "allele frequency", y = "SNP count") +
  theme_light()+
  theme(legend.position = "none",
        plot.subtitle = element_text(hjust = .5))

p3 <- data |> 
  mutate(maf = minor/(minor+major)) |> 
  group_by(ind, maf) |> 
  summarise(n = sum(n)) |> 
  ungroup() |> 
  mutate(rep = map(n, \(x){rep(1, x)})) |> 
  unnest(rep) |> 
  ggplot(aes(x = maf)) +
  geom_histogram(binwidth = .025, boundary = 0,
                 color = "gray70", fill = "gray95") +
  geom_vline(data = tibble(a = c(.5, .3,.2)),
             aes(xintercept =  a, linetype = factor(a)),
             alpha = .5, linewidth = .4) +
  scale_linetype_manual(values = c(`0.5` = 1, `0.3` = 2, `0.2` = 3),
                        guide = "none") +
  labs(subtitle = "2D allele frequency spectrum", x = "major allele count", y = "minor allele count") +
  facet_wrap(ind ~ ., scales = "free_y") +
  theme_light()


p_out <- (p3 / p1 + plot_layout(heights = c(1, .15)) ) | p2

scl <- 1.4
ggsave(plot = p_out,
      filename = out_plt,
      width = 15*scl,
      height = 6*scl,
      device = cairo_pdf)