library(tidyverse)
library(here)
library(glue)
library(prismatic)
library(ggridges)
library(patchwork)
source(here("code/R/project_defaults.R"))

on_x <- read_tsv(here("results/genomes/sex_chrom/mirang_sex_chrom.bed"))

pl_h <- 0
pl_wh <- 2
pl_n <- 10
pl_wn <- 50
pl_l <- 10
pl_g <- 1


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


genomes <- specs |> map_dfr(import_genomes)
g_starts <- dplyr::select(genomes, chr, ref = spec, start_pos, eo) 

plink_data <- read_table(here(glue("results/roh/plink/mirang_filtered_all_h{pl_h}_wh{pl_wh}_n{pl_n}_wn{pl_wn}_l{pl_l}_g{pl_g}/mirang_filtered_all_h{pl_h}_wh{pl_wh}_n{pl_n}_wn{pl_wn}_l{pl_l}_g{pl_g}.hom"))) |> 
  mutate(CHR = str_remove(SNP1, ":.*"),
         autosome = !(CHR %in% on_x$chr)) |> 
  left_join(g_starts |> filter(ref == "mirang"), by = c(CHR = "chr")) |> 
  left_join(samples, by = c(IID = "sample")) |> 
  mutate(gstart = POS1 + start_pos,
         gend = POS2 + start_pos,
         gmid = (gstart + gend) / 2 )

p1 <- plink_data |> 
  ggplot(aes(color = autosome)) +
  geom_rect(data = genomes |> filter(spec == "mirang"), 
            aes(xmin = start_pos,
                xmax = end_pos, 
                ymin = -Inf,
                ymax = Inf, 
                fill = as.character(eo)),
            color = "transparent") +
  geom_point(aes(x = gmid, y = IID, fill = after_scale(clr_lighten(color))), size = .75, shape = 21) +
  geom_linerange(aes(xmin = gstart, xmax = gend, y = IID), lwd = 15) +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') + 
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank())
  
p2 <- plink_data |> 
  ggplot(aes(y = IID, x = log10(KB), color = spec)) +
  geom_density_ridges2(aes(fill = after_scale(clr_alpha(color)),
                           group = IID), scale = .7)+ 
  scale_color_manual(values = clrs) + 
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank())

sum_bp <- genomes |> 
  filter(spec == "mirang") |> 
  pluck("length") |> 
  sum()

roh_summary_by_sample <- plink_data |>
  group_by(IID, spec) |> 
  summarize(sum_roh_length = sum(KB * 1e3)) |> 
  # left_join(mappable_lengths |> filter(roh_type == roh_version)) |> 
  ungroup() |> 
  mutate(f_roh = sum_roh_length / sum_bp) # |>
  # filter(Length_bp > roh_length_threshold) |> 

p3 <- roh_summary_by_sample |> 
  ggplot(aes(x = spec, color = spec)) +
  ggbeeswarm::geom_beeswarm(aes(y = f_roh), alpha = .7) +
  scale_color_manual(values = clrs) + 
  theme_minimal(base_family = fnt_sel, base_size = 9)

p1 + p2 + p3 + 
  plot_layout(widths = c(1, .3, .3), nrow = 1, guides = "collect") +
  plot_annotation(title = glue("ROH sumary (max certain, plink)"),
                  caption = glue("only ROH on autosomes are included")) & 
  scale_color_manual(values = c(clrs, `TRUE` = "darkgray", `FALSE` = "red")) &
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        strip.text.y.left = element_text(angle = 0),
        plot.caption = element_text(),
        plot.subtitle = element_text())

ggsave(here(glue("results/img/roh/mirang_roh_max_plink_h{pl_h}_wh{pl_wh}_n{pl_n}_wn{pl_wn}_l{pl_l}_g{pl_g}.pdf")),
       width = 16, height = 6, device = cairo_pdf)
