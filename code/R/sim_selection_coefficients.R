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
  geom_vline(xintercept = c(c(.001, .33)),
             linetype = 3, color = "red") +
  geom_text(data = tibble(x = 10^c(-3.2, -1.75, -0.25),
                          lab = s_clases),
            aes(x = x, y = 4, label = lab, group = lab),
            color = "red", family = fnt_sel) +
  geom_density(color = clr_alpha(clr_default[[1]], .1),
               fill = "transparent",
               linewidth = .2) +
  scale_x_continuous(trans = scales::log10_trans()) +
  labs(y = "density", subtitle = "overall distribution of selection coefficients") +
  theme_ms()

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
  theme(axis.title.x = element_blank())

p_out <- p1 / p2

ggsave("results/img/load/sim_tally.pdf",
       width = 7, height = 5, 
       device = cairo_pdf)
