library(tidyverse)
library(patchwork)
library(glue)
library(prismatic)
library(here)

source(here("code/R/project_defaults_shared.R"))
specs <-  c("mirang", "mirleo")
spec_names <- c(mirang = "M. angustirostris", mirleo = "M. leonina")

import_genomes <- function(spec){
  read_tsv(here("data", "genomes", "filtered", glue("{spec}_filt.fa.gz.fai")),
           col_names = c("chr", "length", "offset", "linebases", "linewidth")) |> 
    mutate(spec = spec,
           idx = row_number(),
           end_pos = cumsum(length),
           start_pos = lag(end_pos, default = 0),
           eo = idx %% 2)
}

read_pop <- \(spec){read_tsv(here(glue("results/pop/inds_{spec}.pop")),col_types = "c", col_names = "sample") |> mutate(spec = spec)}
pops <- c("mirang", "mirleo") |> map_dfr(read_pop) 

genomes <- specs |> map_dfr(import_genomes)

g_starts <- dplyr::select(genomes, chr, ref = spec, start_pos) 

read_load <- \(sample, type){
  read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{type}/{sample}_{type}.bed.gz")),
           col_types = "ciic") |> 
    left_join(g_starts, by = c(`#CHROM` = "chr")) |> 
    mutate(gpos = FROM + start_pos) |> 
    select(chr = `#CHROM`, pos = FROM, gpos) |> 
    mutate(sample = sample, type = type)
}

samples <- dir("results/mutation_load/snp_eff/by_ind/expressed/") |> str_remove("_expressed.bed.gz")

data <- crossing(sample = samples,
                 type = c("expressed_anc", "masked", "fixed_anc")) |> 
  pmap_dfr(read_load) |> 
  left_join(pops)


load_labs <- c(masked = "inbreeding\nload",
               expressed = "segregating", 
               fixed = "drift")

data_counts <- data |> 
  group_by(chr, pos, gpos, type, spec) |> 
  count() |> 
  ungroup() |> 
  mutate(type = str_remove(type, "_anc"),
         load_label = factor(load_labs[type], levels = load_labs[c(2,1,3)]))

genome_mirang <- genomes |> 
  filter(spec == "mirang") |>
  mutate(mid_pos = (start_pos + end_pos) / 2)

clr_load_lab <- clr_load |> 
  set_names(nm = c(total = "total",
                   masked = load_labs[[1]],
                   fixed = load_labs[[3]],
                   expressed = load_labs[[2]])[names(clr_load)])

ggplot() +
  geom_rect(data = genome_mirang , 
            aes(xmin = start_pos,
                xmax = end_pos, 
                ymin = -Inf,
                ymax = Inf, 
                fill = as.character(eo)),
            color = "transparent") +
  geom_point(data = data_counts,
             mapping =  aes(x = gpos,
                            y = n*2*(1.5- as.numeric(factor(spec))),
                            #color = load_label
                            color=  spec
                            ),
             size = .6, alpha = .5) +
  geom_hline(yintercept = 0, linewidth = .3, linetype = 3) +
  scale_x_continuous(name = "Genomic Position (Gb)",
                     labels = \(x){sprintf("%.1f", x * 1e-9)},
                     sec.axis = sec_axis(name = "Scaffold Id",
                                         trans = identity,
                                         breaks = genome_mirang$mid_pos[1:17],
                                         label = genome_mirang$chr[1:17] |>
                                           str_remove("N[CW]_0723")))+
  scale_y_continuous("Sample count at SNP (n)",
                     labels = \(y){abs(y)}#,
                     # sec.axis = sec_axis(name = "Load Frequency",
                     #                     trans = identity,
                     #                     breaks = c(-20,-15,-10,-5,0,5,10,15,20),
                     #                     label = c(1,0.75, 0.5, 0.25, 0,
                     #                               0.25, 0.5, 0.75, 1))
                     ) +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
  # scale_color_manual(values = clr_load_lab) +
  scale_color_manual("Species",
                     values = clr_default,
                     labels = spec_names) +
  facet_grid(load_label ~ ., switch = "y") +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  coord_cartesian(ylim = c(-21, 21), expand = 0) +
  theme_ms() +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.text = element_text(face = "italic"))

ggsave("results/img/load_along_genome_draft.pdf", width = 12, height = 5, device = cairo_pdf)
