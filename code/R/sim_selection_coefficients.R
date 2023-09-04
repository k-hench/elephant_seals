library(tidyverse)
library(glue)
library(here)
library(prismatic)
library(patchwork)
library(ggdist)
source(here("code/R/project_defaults_shared.R"))

read_s <- \(idx){
  read_tsv(here(glue("results/simulations/selection_coefficient/sim_{idx}_s.tsv.gz"))) |> 
    mutate(sim_idx = idx)
}

read_tally <- \(idx){
  read_tsv(here(glue("results/simulations/tally/sim_{idx}_load_by_ind.tsv"))) |> 
    mutate(sim_idx = idx)
}

s_clases <- c("weak", "intermediate", "strong")

data <- map_dfr(1:100, read_s )
data_tally <- map_dfr(1:100, read_tally ) |> 
  pivot_longer(cols = starts_with("n"),
               names_sep = "_", 
               names_to = c("pre", "load_type", "sel_class"),
               values_to = "n")

p1 <- data |> 
  ggplot(aes(x = -S,
             group = sim_idx)) +
  geom_line(data = tibble(x = 10^seq(-3.5,0, length.out = 351),
                          y = dgamma(x, shape = 0.04, rate = 0.2)),
            aes(x = x, y = y *.05, group = "input"),
            linetype = 3, color = clr_default[2])+
  geom_vline(xintercept = c(c(.001, .33)),
             linetype = 3, color = clr_default[1]) +
  geom_text(data = tibble(x = 10^c(-3.3, -1.75, -0.15),
                          lab = s_clases),
            aes(x = x, y = 4.5, label = lab, group = lab),
            color = clr_default[1],
            family = fnt_sel,
            size = 3) +
  geom_density(color = clr_alpha(clr_default[[1]], .1),
               fill = "transparent",
               linewidth = .2) +
  scale_x_continuous(trans = scales::log10_trans()) +
  labs(y = "density", subtitle = "overall distribution of selection coefficients") +
  coord_cartesian(ylim = c(0, 4.65), expand = 0) +
  theme_ms() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

p2 <- data_tally |> 
  mutate(sel_class = factor(sel_class, levels = s_clases)) |> 
  ggplot(aes(x = load_type, y = n))+
  stat_slab(aes(thickness = after_stat(pdf),
                fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, .95, .66)))#,
                # color = load_type,
                # fill = after_scale(color)
                ),
            color = clr_default[[1]],
            fill = clr_default[[1]],
            linewidth = .25, scale = .75, side = "both",
            adjust = .5,
            normalize = "groups",
            trim = FALSE, n = 301) +
  scale_colour_ramp_discrete(from = "white", 
                             aesthetics = "fill_ramp",
                             guide = "none", range = c(1, 0)) +
  stat_dots(side = "both", scale = 0.3, quantiles = 100,
            color = clr_alpha("black"), shape = 19) +
  # scale_color_manual(values = clr_load,
  #                    guide = "none") +
  labs(subtitle = "individual load scores (n SNPs)") +
  facet_wrap(sel_class ~ ., scales = "free")+
  theme_ms() +
  theme(axis.title.x = element_blank(),
        strip.text = element_text(color = clr_default[1]))

p_out <- p1 / p2 +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = fnt_sel))

ggsave("results/img/load/sim_tally.pdf",
       width = 7, height = 5, 
       device = cairo_pdf)
