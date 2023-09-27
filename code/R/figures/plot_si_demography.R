library(tidyverse)
library(ggdist)
library(ggtext)
library(here)
library(glue)
library(patchwork)
source(here("code/R/project_defaults_shared.R"))

demtype_ord <- c("whg_null","whg_bot10" ,"whg_bot06", "rad_null",  "rad_bot10", "rad_bot06", "rad_bot06_1k")

data_boot <- read_tsv(here("results/demography/all_models_raw_boots.tsv")) |>  mutate(demtype = factor(demtype, levels = demtype_ord))
data_ci <- read_tsv(here("results/demography/all_models_ci_boots.tsv")) |>  mutate(demtype = factor(demtype, levels = demtype_ord))
data_estimates <- read_tsv(here("results/demography/all_models_estimates.tsv")) |>  mutate(demtype = factor(demtype, levels = demtype_ord))

target_dem <- "rad_bot06_1k"

plot_param <- \(param){
ggplot(mapping = aes(x = .data[[param]], y = demtype)) +
  stat_slab(data = data_boot,
            aes(color = demtype == target_dem,
                fill = after_scale(clr_alpha(color))),
            normalize = "xy",
            adjust = .8,
            linewidth = .2,
            height = .8,
            trim = FALSE, density = density_unbounded ) +
  geom_linerange(data = data_ci |> filter(stat == param),
                 aes(x = NULL, xmin = ci_66_l, xmax = ci_66_u),
                 linewidth = 1.5) +
  geom_linerange(data = data_ci |> filter(stat == param),
                 aes(x = NULL, xmin = ci_95_l,xmax = ci_95_u)) +
  geom_point(data = data_estimates,
             shape = 21,
             size = 3,
             stroke = .4,
             aes(color = demtype == target_dem), 
             fill = clr_alpha("white")) +
  scale_color_manual(values = c(`TRUE` = clr_darken(clr_default[[1]]),
                                `FALSE` = clr_darken(clr_default[[2]])),
                     guide = "none") +
   labs(x = p_units[[param]]) +
  theme_ms() +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown())
}

p_units <- c(NLGM = "*N<sub>eLGM</sub>* (*N<sub>e</sub>*)",
             NANC = "*N<sub>ePREBOT</sub>* (*N<sub>e</sub>*)",
             NBOT = "*N<sub>eBOT</sub>* (*N<sub>e</sub>*)",
             NCUR = "*N<sub>ePOSTBOT</sub>* (*N<sub>e</sub>*)",
             T1 = "*T<sub>se</sub>* (generations)")

p_out <- names(p_units) |> 
  map(plot_param) |> 
  wrap_plots() +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")") &
  theme(text = element_text(family = fnt_sel))

ggsave(filename = here("results/img/final/sf_dem.pdf"),
       plot = p_out,
       width = 8, height = 5.5,
       device = cairo_pdf)