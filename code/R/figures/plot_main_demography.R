library(tidyverse)
library(ggdist)
library(ggtext)
library(here)
library(glue)
library(patchwork)
source(here("code/R/project_defaults_shared.R"))

# --- set key variable ---------
target_dem <- "rad_bot06_1k"
ne_levels <- c("NCUR", "NBOT", "NANC", "NLGM")
ne_labels <- c(NLGM = "*N*<sub>eLGM</sub>",
               NANC = "*N*<sub>ePREBOT</sub>",
               NBOT = "*N*<sub>eBOT</sub>",
               NCUR = "*N*<sub>ePOSTBOT</sub>")

clr_range <- "black"

# --- helper functions -------
geom_ci <- \(data_e = data_estimates,
             data_c = data_ci,
             color = clr_range,
             lwd = c(.5, 1.25)){
  list(
    geom_linerange(data = data_c,
                   inherit.aes = FALSE,
                   aes(xmin = ci_95_l,
                       xmax = ci_95_u,
                       y = stat),
                   linewidth = lwd[1],
                   color = color),
    geom_linerange(data = data_c,
                   inherit.aes = FALSE,
                   aes(xmin = ci_66_l,
                       xmax = ci_66_u,
                       y = stat),
                   linewidth = lwd[2],
                   color = color),
    geom_point(data = data_e,
               shape = 21,
               size = 3,
               stroke = lwd[1],
               color = color, 
               fill = "white")
  )
}


# --- read in data -----------
data_boot <- read_tsv(here("results/demography/all_models_raw_boots.tsv")) |> 
  filter(demtype == target_dem) |> 
  select(starts_with("N")) |> 
  pivot_longer(everything(), names_to = "stat") |> 
  mutate(stat = factor(stat, levels = ne_levels))

data_ci <- read_tsv(here("results/demography/all_models_ci_boots.tsv")) |>
  filter(demtype == target_dem) |> 
  mutate(stat = factor(stat, levels = ne_levels))

data_estimates <- read_tsv(here("results/demography/all_models_estimates.tsv")) |>
  filter(demtype == target_dem) |> 
  select(starts_with("N")) |> 
  pivot_longer(everything(), names_to = "stat") |> 
  mutate(stat = factor(stat, levels = ne_levels))

p1 <- data_boot |> 
  ggplot(aes(x = value, y = stat)) +
  geom_path(data = data_estimates |> 
              arrange(stat),
            aes(group = 1),
            linetype = 2,
            linewidth = .3,
            color = clr_red_line) +
  stat_slab(data = data_boot,
            color = clr_default[[1]],
            fill = clr_alpha(clr_default[[1]]),
            normalize = "xy",
            adjust = .8,
            linewidth = .5,
            height = .9,
            trim = FALSE, density = density_unbounded ) +
  geom_ci() +
  scale_y_discrete(limits = factor(ne_levels), label = \(x){glue("{ne_labels[x]}")}) +
  scale_x_continuous(breaks = (0:5)*5000) +
  coord_cartesian(xlim = c(0, 28000)) +
  labs(x = "*N*<sub>e</sub>", y = "Demographic parameter") +
  theme_ms() +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        axis.title.x = ggtext::element_markdown(family = fnt_sel),
        axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.text.y = ggtext::element_markdown(family = fnt_sel),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(hjust = .5))

p2 <- data_boot |> 
  filter(stat %in% c("NBOT")) |> 
  ggplot(mapping = aes(y = stat, x = value)) +
  stat_slab(color = clr_default[[1]],
            fill = clr_alpha(clr_default[[1]]),
            adjust = 2.5, height = .6,
            linewidth = .5,
            trim = FALSE, density = density_unbounded ) +
  geom_ci(data_e = data_estimates |> 
            filter(stat %in% c("NBOT")),
          data_c = data_ci |> 
            filter(stat %in% c("NBOT"))) +
  labs(x = "*N<sub>e</sub>*", x = NULL) +
  coord_cartesian(ylim = c(.45, 1.65),
             xlim = c(4, 8.5),
             expand = 0) +
  theme_ms() +
  theme(plot.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        # axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(hjust = .5))

p3 <- data_boot |> 
  filter(stat %in% c("NCUR")) |> 
  ggplot(mapping = aes(x = value, y = stat)) +
  stat_slab(color = clr_default[[1]],
            fill = clr_alpha(clr_default[[1]]),
            adjust = 2.5, height = .6,
            linewidth = .5,
            trim = FALSE, density = density_unbounded ) +
  geom_ci(data_e = data_estimates |> 
            filter(stat %in% c("NCUR")),
          data_c = data_ci |> 
            filter(stat %in% c("NCUR"))) +
  labs(x = "*N<sub>e</sub>*", y = NULL) +
  coord_cartesian(ylim = c(.45, 1.55),
                  xlim = c(2455, 2855),
                  expand = 0) +
  theme_ms() +
  theme(plot.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        # axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(hjust = .5))


p_out <- p1 +
    annotation_custom(grob = ggplotGrob(p2),
                      ymin = 1.475,
                      ymax = 2.4,
                      xmin = 14500, 
                      xmax = 27000) +
    annotation_custom(grob = ggplotGrob(p3),
                      ymin = 0.425,
                      ymax = 1.4,
                      xmin = 14500, 
                      xmax = 27000)

ggsave(filename = here("results/img/final/f_dem.pdf"),
       plot = p_out,
       width = 8, height = 5.5,
       device = cairo_pdf)

ggsave(filename = here("results/img/final/f_dem.png"),
       plot = p_out,
       width = 8, height = 5.5)
