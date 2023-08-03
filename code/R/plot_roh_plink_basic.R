library(tidyverse)
library(here)
library(glue)
library(prismatic)
library(ggridges)
library(patchwork)
source(here("code/R/project_defaults.R"))

on_x <- read_tsv(here("results/genomes/sex_chrom/mirang_sex_chrom.bed"))

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

plink_path <- here("results/roh/plink/")
plink_folders <- dir(plink_path, pattern = "d[0-9]*$")

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
  
  p_out <- p1 + p2 + p3 + 
    plot_layout(widths = c(1, .3, .3), nrow = 1, guides = "collect") +
    plot_annotation(title = glue("ROH sumary (max certain, plink)"),
                    subtitle = glue("window-snp {pl_wn} snp {pl_n} kb {pl_l} gap {pl_g} density {pl_den} window-missing {pl_mis} het {pl_h} window-het {pl_wh}"),
                    caption = glue("only ROH on autosomes are included")) & 
    scale_color_manual(values = c(clrs, `TRUE` = "darkgray", `FALSE` = "red")) &
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          strip.text.y.left = element_text(angle = 0),
          plot.caption = element_text(),
          plot.subtitle = element_text(family = "ubuntu mono"))
  
  ggsave(here(glue("results/img/roh/mirang_roh_max_plink_{plink_tag}.pdf")),
         plot = p_out,
         width = 16, height = 6, device = cairo_pdf)
}

plink_folders |> walk(summarize_plink_run)
