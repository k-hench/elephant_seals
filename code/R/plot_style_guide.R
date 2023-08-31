library(tidyverse)
library(prismatic)
library(patchwork)
library(here)
library(ggdist)
library(extrafont)
library(cowplot)
#font_import()
loadfonts()
source(here("code/R/project_defaults_shared.R"))

p_traj <- tibble(x = 0:51) |> 
  ggplot() +
  geom_tile(aes(x = x, y = 1, fill = x)) +
  coord_fixed(ratio = 1.5) +
  scale_fill_gradientn(colours = clr_f_traj(7),
                       guide = "none")+
  labs(subtitle = "color scheme trajectories {clr_f_traj()}") +
  theme_void(base_family = fnt_sel) +
  theme(axis.text.x = element_text(),
        plot.subtitle = element_text(hjust = .5))

p_pheno <- tibble(pheno = names(clr_pheno)) |> 
  ggplot() +
  geom_point(aes(x = pheno, y = .5,
                 fill = pheno, color = after_scale(clr_darken(fill, .6))),
             size = 10, shape = 21) +
  coord_fixed(ratio = 1,
              xlim = c(.5, length(names(clr_pheno))+.5),
              ylim = c(.25,.75)) +
  scale_fill_manual(values = clr_pheno, guide = "none") +
  labs(subtitle = "color scheme phenotypes {clr_pheno}") +
  theme_void(base_family = fnt_sel) +
  theme(axis.text.x = element_text(),
        plot.subtitle = element_text(hjust = .5))

p_load <- tibble(load = factor(names(clr_load),
                               levels = names(clr_load))) |> 
  ggplot() +
  geom_point(aes(x = load, y = .5,
                 fill = load, color = after_scale(clr_darken(fill, .6))),
             size = 10, shape = 21) +
  coord_fixed(ratio = 1,
              xlim = c(.5, length(names(clr_load))+.5),
              ylim = c(.25,.75)) +
  scale_fill_manual(values = clr_load, guide = "none") +
  labs(subtitle = "color scheme  loadtypes {clr_load}") +
  theme_void(base_family = fnt_sel) +
  theme(axis.text.x = element_text(),
        plot.subtitle = element_text(hjust = .5))

p_def <- tibble(type = c("default", "secondary")) |> 
  ggplot() +
  geom_point(aes(x = type, y = .5,
                 fill = type, color = after_scale(clr_darken(fill, .6))),
             size = 10, shape = 21) +
  coord_fixed(ratio = 1,
              xlim = c(.5, 2.5),
              ylim = c(.25,.75)) +
  scale_fill_manual(values = clr_default |> set_names(nm = c("default", "secondary")),
                    guide = "none") +
  labs(subtitle = "default color scheme {clr_default}") +
  theme_void(base_family = fnt_sel) +
  theme(axis.text.x = element_text(),
        plot.subtitle = element_text(hjust = .5))

p_scatter <- ggplot(iris, aes(x = Petal.Length, y = Sepal.Length)) +
  geom_point(aes(color = after_scale(clr_darken(fill, .6)),
                   fill = Species == "versicolor"),
               shape = 21)  +
  scale_fill_manual(values = c(`TRUE` = clr_default[[1]],  `FALSE` = clr_default[[2]]),
                    guide = "none") +
  facet_wrap(Species ~ . ) +
  theme_ms()

p_b <- ggplot(iris, aes(x = Species, y = Sepal.Length)) +
  geom_boxplot(aes(color = after_scale(clr_darken(fill, .6)),
               fill = Species == "versicolor"))  +
  scale_fill_manual(values = c(`TRUE` = clr_default[[1]],  `FALSE` = clr_default[[2]]),
                    guide = "none") +
  theme_ms()

p_v <- ggplot(iris, aes(x = Species, y = Sepal.Length)) +
  stat_slab(aes(thickness = after_stat(pdf),
                fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, .95, .66))),
                color = Species == "versicolor",
                fill = after_scale(color)),
            linewidth = .25, scale = .75, side = "both",
            adjust = .5,
            normalize = "groups",
            trim = FALSE, n = 301) +
  scale_colour_ramp_discrete(from = "white", 
                             aesthetics = "fill_ramp",
                             guide = "none", range = c(1, 0)) +
  stat_dots(side = "both", scale = 0.3, quantiles = 100,
            color = clr_alpha("black"), shape = 19) +
  scale_color_manual(values = c(`TRUE` = clr_default[[1]],
                                `FALSE` = clr_default[[2]]),
                    guide = "none") +
  theme_ms()

p_gg <- ggplot(iris,
       aes(x = Sepal.Length, y = Sepal.Width, group = Species)) +
  geom_smooth(method = "lm", se = FALSE, aes(color = Species == "versicolor"),
              linewidth = .5, fullrange = TRUE) +
  geom_point(aes(color = after_scale(clr_darken(fill, .6)),
                 fill = Species == "versicolor"),
             shape = 21)  +
  scale_color_manual(values = c(`TRUE` = clr_default[[1]],  `FALSE` = clr_default[[2]]),
                    guide = "none")  +
  scale_fill_manual(values = c(`TRUE` = clr_default[[1]],  `FALSE` = clr_default[[2]]),
                    guide = "none")  +
  coord_cartesian(xlim = range(iris$Sepal.Length),
                  ylim = range(iris$Sepal.Width)) +
  labs(subtitle = "ggplot") +
  theme_ms()

p1 <- function() {
  par(
    mar = c(3, 3, 1, 1),
    mgp = c(2, 1, 0)
  )
  plot(Sepal.Width ~ Sepal.Length, 
       pch = 21,
       bg = clr_default[2 - (iris$Species == "versicolor")],
       col = clr_darken(clr_default[2 - (iris$Species == "versicolor")]),
       data = iris,
       family = fnt_sel,
       bty = "l",
       tcl = -.33,
       cex =.6,
       main = "base R plot",
       family = fnt_sel,
       cex.axis = 0.75, las = 1, tcl = -0.25)
  abline(lm(Sepal.Width ~ Sepal.Length, data = iris[iris$Species == "versicolor",]),
         col = clr_default[1])
  abline(lm(Sepal.Width ~ Sepal.Length, data = iris[iris$Species == "virginica",]),
         col = clr_default[2])
  abline(lm(Sepal.Width ~ Sepal.Length, data = iris[iris$Species == "setosa",]),
         col = clr_default[2])
}

plot_grid(p_def,p_load, p_pheno, p_traj, ncol = 1) +
plot_grid(p_scatter) /
plot_grid(p_b, p_v) /
plot_grid(p1, p_gg) +
  plot_layout(widths = c(.7, 1))

scl <- 1.1
ggsave("~/Dropbox/Elephant seal paper/figures/drafts_and_discussion/style_guide.pdf", 
       width = 16 * scl,
       height = 9 * scl,
       device = cairo_pdf)
