library(tidyverse)
library(prismatic)
library(here)
library(glue)
library(ggdist)
source(here("code/R/project_defaults_shared.R"))

q_quant <- c(.95, 0.66)
q_borders <- c((1 - q_quant[1])/2, (1 - q_quant[2])/2)

rad_tags <- c("Bott10", "Bott6_Best", "Null")
rad_translate <- str_c("RAD_",
                       c("bot10", "bot06", "null")) |> 
  str_c("-lgm") |> 
  set_names(rad_tags)

estimates_rad <- tribble(
  ~stat,	~Estimate,	~lowerCI, ~upperCI,
  "NGM",	266.5,	226.5,	1721.5,
  "T1",	1222,	350,	1570,
  "NANC",	12855.5,	2827.5,	20274.5,
  "NBOT",	5.5,	5,	7.5,
  "NCUR",	2624,	2506.5,	2773 ) |> 
  mutate(demtype = "RAD_bot06-lgm")

read_rad <- \(tag){
  read_tsv(here(glue("results/RAD/demography/Non_Param_boots_{tag}Model.txt"))) |> 
    mutate(demtype = rad_translate[tag])
}

data_boot_rad <- read_rad("Bott6_Best") |> 
  bind_rows(read_tsv("../david/NonParam_900ExtraBoots.txt") |> 
              mutate(demtype = "RAD_bot06-lgm"))

data_boot_ci <- data_boot_rad |> 
  group_by(demtype) |> 
  summarise(across(everything(), 
                   \(x){
                     list(tibble(mean = mean(x),
                                 median = median(x),
                                 sd = sd(x),
                                 q_in = str_c(quantile(x, c(1 - q_borders[2], q_borders[2]), na.rm = TRUE), collapse = "_"),
                                 q_out = str_c(quantile(x, c(1 - q_borders[1], q_borders[1]), na.rm = TRUE), collapse = "_")))
                   })) |> 
  ungroup() |> 
  pivot_longer(-demtype, names_to = "stat") |> 
  unnest(value) |> 
  separate(q_in, into = c("ci_in_l", "ci_in_u"), sep = "_", convert = TRUE) |> 
  separate(q_out, into = c("ci_out_l", "ci_out_u"), sep = "_", convert = TRUE)

n_levels <- c("NCUR",  "NBOT", "NANC", "NLGM")

mod_labels <- c(NLGM = "*N<sub>eLGM</sub>*",
                NANC = "*N<sub>ePREBOT</sub>*",
                NBOT = "*N<sub>eBOT</sub>*",
                NCUR = "*N<sub>ePOSTBOT</sub>*")

data_est <- estimates_rad |>
  filter(stat != "T1") |> 
  select(demtype, stat, n = Estimate) |>
  mutate(param = factor(str_replace(stat, "NGM", "NLGM"), levels = rev(n_levels)))
  
data_bot_plot <- data_boot_rad |>
  select(demtype, starts_with("N")) |> 
  mutate(demtype = str_replace(demtype, "RAD_", "rad-") |> 
           str_remove("-lgm")) |> 
  rename(NLGM = "NGM") |> 
  pivot_longer(cols = -demtype,
               names_to = "param", 
               values_to = "n",
               names_transform = \(x){factor(x, levels = rev(n_levels))})
  
p0 <- data_bot_plot |> 
  filter(param %in% c("NBOT")) |> 
  ggplot(aes(y = n, x = param)) +
  stat_slab(color = "transparent", 
            fill = "transparent",
            n = 2) +
  geom_linerange(data = data_est |> 
                   filter(param %in% c("NBOT")),
                 linewidth = .3,
                 color = clr_default[[1]],
                 aes(xmin = as.numeric(param)-.45-2,
                     xmax = as.numeric(param)+.45-2)) +
  stat_slab(aes(thickness = after_stat(pdf),
                fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, q_quant)))),
            color = clr_default[[1]],
            fill  = clr_default[[1]],
            linewidth = .45, scale = .75, side = "both",
            # adjust = .3,
            adjust = 2,
            normalize = "groups",
            trim = FALSE, n = 301) +
  scale_colour_ramp_discrete(from = "white", aesthetics = "fill_ramp", guide = "none", range = c(1,0)) +
  stat_dots(side = "both", scale = 0.3, quantiles = 100,
            color = clr_alpha("black"), shape = 19) +
  geom_point(data = data_est |> 
               filter(param %in% c("NBOT")),
             shape = 21,
             size = 3,
             stroke = .4,
             color = clr_darken(clr_default[[1]]), 
             fill = clr_alpha("white")) +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  coord_cartesian(xlim = c(.45, 1.55),
                  ylim = c(4, 8.5),
                  expand = 0) +
  theme_ms() +
  theme(plot.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.title = element_blank(),
        # axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.subtitle = element_text(hjust = .5))


single_dot_slab <- \(par){
  stat_dots(data = data_bot_plot |> 
              filter(param == par),
            side = "both", scale = 0.4, 
            quantiles = 100,
            color = clr_alpha("black"), 
            shape = 19)
}

p1 <- data_bot_plot |> 
  ggplot(aes(y = n, x = param)) +
  stat_slab(color = "transparent", 
            fill = "transparent",
            n = 2) +
  # geom_vline(xintercept = 3, linetype = 3) +
  # geom_rect(data = tibble(x1 = 2.5,
  #                         x2 = 3.5,
  #                         y1 = -7e2,
  #                         y2 = 7e2),
  #           inherit.aes = FALSE,
  #           aes(xmin = x1, xmax = x2, 
  #               ymin = y1, ymax = y2),
  #           color = "red",
  #           fill= "transparent") 
  geom_linerange(data = data_est,
                 linewidth = .3,
                 color = clr_default[[1]],
                 aes(xmin = as.numeric(param)-.45,
                     xmax = as.numeric(param)+.45)) +
  stat_slab(aes(thickness = after_stat(pdf),
                fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, q_quant)))),
            color = clr_default[[1]],
            fill  = clr_default[[1]],
            linewidth = .45, scale = .75, side = "both",
            adjust = .3,
            normalize = "groups",
            trim = FALSE, n = 301) +
  scale_colour_ramp_discrete(from = "white", aesthetics = "fill_ramp", guide = "none", range = c(1,0)) +
  geom_line(data = data_est,
            aes(x = as.numeric(param)),
            linetype = 3,
            color = clr_darken(clr_default[[1]])) +
  (map(c("NLGM", "NANC", "NBOT", "NCUR"), single_dot_slab)) +
  geom_point(data = data_est,
             shape = 21,
             size = 3,
             stroke = .4,
             color = clr_darken(clr_default[[1]]), 
             fill = clr_alpha("white")) +
  scale_x_discrete(label = \(x){glue("{mod_labels[x]}")}) +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  theme_ms() +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        axis.title.y = ggtext::element_markdown(),
        # axis.title.x = element_blank(),
        axis.text.x = ggtext::element_markdown(family = fnt_sel),
        axis.ticks.x = element_blank(),
        plot.subtitle = element_text(hjust = .5))

p1b <- p1 +
  # geom_segment(data = tibble(x1 = c(2.9, 3.1),
  #                            x2 = c(2.55, 3.45),
  #                            y1 = 7e2,
  #                            y2 = 18e3),
  #              inherit.aes = FALSE,
  #              aes(x = x1, xend = x2, 
  #                  y = y1, yend = y2),
  #              linetype = 3) +
  annotation_custom(grob = ggplotGrob(p0),
                    xmin = 2.275, xmax = 3.65, 
                    ymin = 16000, 
                    ymax = 27000) +
  coord_cartesian(xlim = c(.45, 4.55),
                  ylim = c(-1e3,2.7e4),
                  expand = 0)

ggsave(filename = here("results/img/demography/param_single_panel_nbot_only_b.pdf"),
       plot = p1b, 
       width = 9, height = 6, device = cairo_pdf)

p2 <- data_bot_plot |> 
  filter(param %in% c("NCUR")) |> 
  ggplot(aes(y = n, x = param)) +
  stat_slab(color = "transparent", 
            fill = "transparent",
            n = 2) +
  geom_linerange(data = data_est |> 
                   filter(param %in% c("NCUR")),
                 linewidth = .3,
                 color = clr_default[[1]],
                 aes(xmin = as.numeric(param)-.45-3,
                     xmax = as.numeric(param)+.45-3)) +
  stat_slab(aes(thickness = after_stat(pdf),
                fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, q_quant)))),
            color = clr_default[[1]],
            fill  = clr_default[[1]],
            linewidth = .45, scale = .75, side = "both",
            adjust = .3,
            normalize = "groups",
            trim = FALSE, n = 301) +
  scale_colour_ramp_discrete(from = "white", aesthetics = "fill_ramp", guide = "none", range = c(1,0)) +
  stat_dots(side = "both", scale = 0.3, quantiles = 100,
            color = clr_alpha("black"), shape = 19) +
  geom_point(data = data_est |> 
               filter(param %in% c("NCUR")),
             shape = 21,
             size = 3,
             stroke = .4,
             color = clr_darken(clr_default[[1]]), 
             fill = clr_alpha("white")) +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  coord_cartesian(xlim = c(.45, 1.55),
                  ylim = c(2480, 2850),
                  expand = 0) +
  theme_ms() +
  theme(plot.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.title = element_blank(),
        # axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.subtitle = element_text(hjust = .5))

p1c <- p1 +
  # geom_segment(data = tibble(x1 = c(3, 4),
  #                            x2 = c(3, 4),
  #                            y1 = c(7e2, 3.5e3),
  #                            y2 = c(16.5e3, 8.7e3)),
  #              inherit.aes = FALSE,
  #              aes(x = x1, xend = x2, 
  #                  y = y1, yend = y2),
  #              linetype = 3) +
  annotation_custom(grob = ggplotGrob(p0),
                    xmin = 2.275 + .1,
                    xmax = 3.65 - .13, 
                    ymin = 16000, 
                    ymax = 25500) +
  annotation_custom(grob = ggplotGrob(p2),
                    xmin = 3.25, xmax = 4.565, 
                    ymin = 8300, 
                    ymax = 17300) +
  coord_cartesian(xlim = c(.45, 4.55),
                  ylim = c(-1e3,2.7e4),
                  expand = 0)

ggsave(filename = here("results/img/demography/param_single_panel_nbot_ncur_b_update.pdf"),
       plot = p1c,
       width = 9, height = 6, device = cairo_pdf)



# ===== Beeswarm =======
geom_ci <- \(data = data_est,
             data_ci = data_boot_ci_n,
             color = clr_darken(clr_range[[1]]),
             lwd = c(.5, 1.25)){
  list(
    geom_linerange(data = data_ci,
                   inherit.aes = FALSE,
                   aes(ymin = ci_out_l,
                       ymax = ci_out_u,
                       x = stat),
                   linewidth = lwd[1],
                   color = color),
      geom_linerange(data = data_ci,
                     inherit.aes = FALSE,
                     aes(ymin = ci_in_l,
                         ymax = ci_in_u,
                         x = stat),
                     linewidth = lwd[2],
                     color = color),
      geom_point(data = data,
                 shape = 21,
                 size = 3,
                 stroke = lwd[1],
                 color = color, 
                 fill = "white")
  )
}

data_boot_ci_n <- data_boot_ci |> 
  filter(grepl("^N", stat)) |> 
  mutate(stat = if_else(stat == "NGM", "NLGM", stat))

j_size <- .6
p3a <- data_bot_plot |> 
  ggplot(aes(y = n, x = param)) +
  geom_jitter(width = .25, height = 0, size = j_size,
              color = clr_alpha(clr_default[[1]], .3)) +
  # geom_swarm(color = "transparent") +
  # (map(c("NLGM", "NANC", "NBOT", "NCUR"), single_dot_slab)) +
  geom_line(data = data_est,
            aes(x = 5 - as.numeric(param)),
            linetype = 3,
            color = clr_darken(clr_default[[1]])) +
  geom_ci(color = clr_darken(clr_default[[1]])) +
  scale_x_discrete(limits = factor(n_levels), label = \(x){glue("{mod_labels[x]}")}) +
  coord_flip() +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  theme_ms() +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.title.x = ggtext::element_markdown(family = fnt_sel),
        axis.text.y = ggtext::element_markdown(family = fnt_sel),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.subtitle = element_text(hjust = .5))

p3b <- data_bot_plot |> 
  filter(param %in% c("NBOT")) |> 
  ggplot(mapping = aes(y = n, x = param)) +
  geom_jitter(width = .25, height = 0, size = j_size,
              color = clr_alpha(clr_default[[1]], .3)) +
    geom_ci(data = data_est |> 
            filter(param %in% c("NBOT")),
          data_ci = data_boot_ci_n |> 
            filter(stat %in% c("NBOT")),
          color = clr_darken(clr_default[[1]])) +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  coord_flip(xlim = c(.45, 1.55),
                  ylim = c(4, 8.5),
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

p3c <- data_bot_plot |> 
  filter(param %in% c("NCUR")) |> 
  ggplot(mapping = aes(y = n, x = param)) +
  geom_jitter(width = .25, height = 0, size = j_size,
              color = clr_alpha(clr_default[[1]], .3)) +
  geom_ci(data = data_est |> 
            filter(param %in% c("NCUR")),
          data_ci = data_boot_ci_n |> 
            filter(stat %in% c("NCUR")),
          color = clr_darken(clr_default[[1]])) +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  coord_flip(xlim = c(.45, 1.55),
                  ylim = c(2480, 2850),
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


(p3 <- p3a +
  # geom_segment(data = tibble(x1 = c(3, 4),
  #                            x2 = c(3, 4),
  #                            y1 = c(7e2, 3.5e3),
  #                            y2 = c(16.5e3, 8.7e3)),
  #              inherit.aes = FALSE,
  #              aes(x = x1, xend = x2, 
  #                  y = y1, yend = y2),
  #              linetype = 3) +
  # geom_vline(linetype = 2, xintercept = c(1,2)) +
  annotation_custom(grob = ggplotGrob(p3b),
                    xmin = 1.63,
                    xmax = 2.25,
                    ymin = 16000, 
                    ymax = 25500) +
  annotation_custom(grob = ggplotGrob(p3c),
                    xmin = 0.63,
                    xmax = 1.25,
                    ymin = 16000, 
                    ymax = 25500)
)

ggsave(filename = here("results/img/demography/param_single_panel_nbot_ncur_b_alternative.pdf"),
       plot = p3,
       width = 8, height = 4,
       device = cairo_pdf)
# ======================
clr_range <- clr_default[[1]]

p4a <- data_bot_plot |> 
  ggplot(aes(y = n, x = param)) +
  # geom_swarm(color = "transparent") +
  # (map(c("NLGM", "NANC", "NBOT", "NCUR"), single_dot_slab)) +
  geom_line(data = data_est,
            aes(x = 5 - as.numeric(param)),
            linetype = 3,
            color = clr_range) +
  geom_ci(color = clr_range) +
  scale_x_discrete(limits = factor(n_levels), label = \(x){glue("{mod_labels[x]}")}) +
  coord_flip() +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  theme_ms() +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.title.x = ggtext::element_markdown(family = fnt_sel),
        axis.text.y = ggtext::element_markdown(family = fnt_sel),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.subtitle = element_text(hjust = .5))

p4b <- data_bot_plot |> 
  filter(param %in% c("NBOT")) |> 
  ggplot(mapping = aes(y = n, x = param)) +
  geom_ci(data = data_est |> 
            filter(param %in% c("NBOT")),
          data_ci = data_boot_ci_n |> 
            filter(stat %in% c("NBOT")),
          color = clr_range) +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  coord_flip(xlim = c(.45, 1.55),
             ylim = c(4, 8.5),
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

p4c <- data_bot_plot |> 
  filter(param %in% c("NCUR")) |> 
  ggplot(mapping = aes(y = n, x = param)) +
  geom_ci(data = data_est |> 
            filter(param %in% c("NCUR")),
          data_ci = data_boot_ci_n |> 
            filter(stat %in% c("NCUR")),
          color = clr_range) +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  coord_flip(xlim = c(.45, 1.55),
             ylim = c(2480, 2850),
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


p4 <- p4a +
    # geom_vline(linetype = 2, xintercept = c(1,2), linewidth = .2) +
    annotation_custom(grob = ggplotGrob(p4b),
                      xmin = 1.5,
                      xmax = 2.25,
                      ymin = 10000, 
                      ymax = 17500) +
    annotation_custom(grob = ggplotGrob(p4c),
                      xmin = 0.5,
                      xmax = 1.25,
                      ymin = 10000, 
                      ymax = 17500)

ggsave(filename = here("results/img/demography/param_single_panel_nbot_ncur_b_alternative2.pdf"),
       plot = p4,
       width = 8, height = 3.5,
       device = cairo_pdf)


# ======================================

p5a <- data_bot_plot |> 
  ggplot(aes(y = n, x = param)) +
  (map(c("NLGM", "NANC", "NBOT", "NCUR"),
       \(x){
         stat_slab(data = data_bot_plot |> filter(param == x),
                   fill = clr_alpha(clr_default[[1]]),
                   trim = FALSE, density = density_unbounded)
         })) +
  geom_line(data = data_est,
            aes(x = 5 - as.numeric(param)),
            linetype = 3,
            color = clr_range) +
  geom_ci(color = clr_darken(clr_default[[1]])) +
  scale_x_discrete(limits = factor(n_levels), label = \(x){glue("{mod_labels[x]}")}) +
  coord_flip() +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  theme_ms() +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.title.x = ggtext::element_markdown(family = fnt_sel),
        axis.text.y = ggtext::element_markdown(family = fnt_sel),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.subtitle = element_text(hjust = .5))

p5b <- data_bot_plot |> 
  filter(param %in% c("NBOT")) |> 
  ggplot(mapping = aes(y = n, x = param)) +
  stat_slab(fill = clr_alpha(clr_default[[1]]),
            adjust = 2.5, width = .66,
            trim = FALSE, density = density_unbounded )+
  geom_ci(data = data_est |> 
            filter(param %in% c("NBOT")),
          data_ci = data_boot_ci_n |> 
            filter(stat %in% c("NBOT")),
          color = clr_darken(clr_default[[1]])) +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  coord_flip(xlim = c(.45, 1.65),
             ylim = c(4, 8.5),
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

p5c <- data_bot_plot |> 
  filter(param %in% c("NCUR")) |> 
  ggplot(mapping = aes(y = n, x = param)) +
  stat_slab(fill = clr_alpha(clr_default[[1]]),
            adjust = 2.5, width = .55,
            trim = FALSE, density = density_unbounded )+
  geom_ci(data = data_est |> 
            filter(param %in% c("NCUR")),
          data_ci = data_boot_ci_n |> 
            filter(stat %in% c("NCUR")),
          color = clr_darken(clr_default[[1]])) +
  labs(y = "*N<sub>e</sub>*", x = NULL) +
  coord_flip(xlim = c(.45, 1.55),
             ylim = c(2455, 2855),
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


(p5 <- p5a +
  # geom_vline(linetype = 2, xintercept = c(1,2), linewidth = .2) +
  annotation_custom(grob = ggplotGrob(p5b),
                    xmin = 1.475,
                    xmax = 2.4,
                    ymin = 14500, 
                    ymax = 27000) +
  annotation_custom(grob = ggplotGrob(p5c),
                    xmin = 0.425,
                    xmax = 1.4,
                    ymin = 14500, 
                    ymax = 27000))

ggsave(filename = here("results/img/demography/param_single_panel_nbot_ncur_b_alternative3.pdf"),
       plot = p5,
       width = 8, height = 5.5,
       device = cairo_pdf)
