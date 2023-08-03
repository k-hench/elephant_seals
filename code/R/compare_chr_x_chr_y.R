library(tidyverse)
library(ggforce)
library(prismatic)
library(here)
library(glue)
library(patchwork)
library(ggtext)
library(plyranges)
library(rlang)
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

get_psl <- \(file){
  vroom::vroom(here::here("results","psl", file),
               delim = "\t",
               col_names = c("matches", "misMatches", "repMatches", "nCount",
                             "qNumInsert", "qBaseInsert", "tNumInsert",
                             "tBaseInsert", "strand", "qName", "qSize", "qStart",
                             "qEnd", "tName", "tSize", "tStart", "tEnd",
                             "blockCount")) |>
    select( tName, tStart, tEnd, qName, qStart, qEnd, strand ) |>
    mutate(tSize = abs(tEnd - tStart),
           qSize = abs(qEnd - qStart))
}


alignment_on_xy <- get_psl("slim_mirang_on_zalcal.psl.gz") |> 
  filter(#tName %in% c("NC_045612.1", "NC_045613.1"),
         qName == "NC_072356.1")
genome <- import_genomes("mirang")

xy_set <- alignment_on_xy |> 
  filter(qSize > 5e3,
         qName == "NC_072356.1") |> 
  left_join(genome, by = c(qName = "chr")) |> 
  mutate(gstart = start_pos + qStart,
         gend = start_pos + qEnd,
         gmid = (gstart + gend)/2) 

ggplot() +
  geom_rect(data = genome |> filter(chr == "NC_072356.1"), 
            aes(xmin = start_pos,
                xmax = end_pos, 
                ymin = -Inf,
                ymax = Inf, 
                fill = as.character(eo)),
            color = "transparent") +
  geom_jitter(data = xy_set,width = 0, height = .1,
              aes(x = gmid, y = tName, color = tName %in% c("NC_045612.1", "NC_045613.1")#tName
                  ),
              size = 1.75,
              pch = 19,
              alpha = .3) +
  geom_linerange(data = xy_set,
                 aes(xmin = gstart, xmax = gend, y = tName, color = tName %in% c("NC_045612.1", "NC_045613.1") #tName
                     ),
                 linewidth = 5) +
  # geom_linerange(aes(xmin = gstart, xmax = gend, y = IID), lwd = 15) +
  coord_cartesian(#ylim = c(.5, 2.5),
                  expand = 0) +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') + 
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank())
