library(tidyverse)
library(prismatic)
library(here)
library(glue)
library(boot)
source(here("code/R/project_defaults.R"))
clr1 <- clr_default[1]
clr2 <- clr_default[2]
clr_accent <- clr1
clr_accent2 <- clr1

basepath <- "results/demography/fastsimcoal/mirang_on_mirang/"

q_quant <- c(.95, 0.66)
q_borders <- c((1 - q_quant[1])/2, (1 - q_quant[2])/2)


read_bs <- \(type, bs_nr){
  bs_idx <- str_pad(bs_nr, width = 2, pad = 0)
  read_tsv(here(glue("{basepath}/{type}/bootstrap/bs_{bs_idx}/bestrun/mirang_on_mirang_{type}.bestlhoods")))
}

bs_stat <- \(type, bs_nr, stat){
  read_bs(type, bs_nr)[stat][[1]]
}

get_parameter_inventory <- \(type){
  names(read_bs(type, 1))
}

read_bootsrap_value <- \(type, stat){
  count_env$n <- (count_env$n + 1L)
  bs_stat(type, count_env$n, stat)
}

read_boot_set <- \(type, stat){
  # reset counter for randomization process in boot::boot()
  count_env$n <- 0L
  current_boot_read <- \(data, i){
    count_env$n <- (count_env$n + 1L)
    data_boots$bs_data[data_boots$demtype == type][[1]][[stat]][[count_env$n]]
  }
  boot_out <- boot(data = "dummy",
                   statistic = current_boot_read, 
                   R = 49)
  
  count_env$n <- 0L
  tibble(demtype = type,
         stat = stat,
         boot = list(boot_out),
         boot_ci = list(boot.ci(boot_out, 
                                conf = q_quant,
                                type = c("norm", "basic"))))
}

get_estimates <- \(type){
  read_tsv(here(glue("{basepath}/{type}/bestrun/mirang_on_mirang_{type}.bestlhoods"))) |> 
    pivot_longer(everything(), names_to = "stat", values_to = "estimate") |> 
    mutate(demtype = type)
}

access_ci <- \(ci, level = 1){str_c(ci$norm[level,][2:3], collapse = "_")}

summarise_bs_data <- \(data){
  data |> 
    select(-bs_dx) |> 
    summarise(across(everything(),
                     \(x){
                       list(tibble(mean = mean(x),
                                   median = median(x),
                                   sd = sd(x),
                                   q90 = str_c(quantile(x, c(1 - q_borders[2], q_borders[2]), na.rm = TRUE), collapse = "_"),
                                   q95 = str_c(quantile(x, c(1 - q_borders[1], q_borders[1]), na.rm = TRUE), collapse = "_")))
                     })) |> 
    ungroup() |> 
    pivot_longer(everything(), names_to = "stat") |> 
    unnest(value) |> 
    separate(q90, into = c("ci_90_l", "ci_90_u"), sep = "_", convert = TRUE) |> 
    separate(q95, into = c("ci_95_l", "ci_95_u"), sep = "_", convert = TRUE)
}

all_dem_types <- c("bot06-lgm", "bot06-nes",
                   "bot10-lgm", "bot10-nes",
                   "null-lgm", "null-nes")

config_table <- tibble(type = all_dem_types) |> 
  mutate(stat = map(type, get_parameter_inventory)) |> 
  unnest(stat)

# read_bs(all_dem_types[1], 4)
# bs_stat(all_dem_types[1], 4, "RTEA")
# get_parameter_inventory(all_dem_types[1])
# set up counter to trick the randomization process in boot::boot()
count_env <- new.env()
count_env$n <- 0L

get_boot_stats <- \(type){
  tibble(demtype = type,
         bs_data = list(0:99 |> 
                          map_dfr(read_bs, type = type) |> 
                          mutate(bs_dx = 0:99))
  )
}

data_boots <- all_dem_types |> 
  map_dfr(get_boot_stats) |> 
  mutate(data_quantiles = map(bs_data, summarise_bs_data))

all_boots <- config_table |> 
  pmap_dfr(read_boot_set) |>
  mutate(ci_90 = map_chr(boot_ci, access_ci, level = 1),
         ci_95 = map_chr(boot_ci, access_ci, level = 2)) |> 
  separate(ci_90, into = c("ci_90_l", "ci_90_u"), sep = "_", convert = TRUE) |> 
  separate(ci_95, into = c("ci_95_l", "ci_95_u"), sep = "_", convert = TRUE)

all_estimates <- all_dem_types |> 
  map_dfr(get_estimates)

# # ci computed by {boot}??
# all_estimates |> 
#   left_join(all_boots |> select(demtype, stat, starts_with("ci_"))) |> 
#   ggplot(aes(y = demtype, color = demtype)) +
#   geom_linerange(aes(xmin = ci_95_l, xmax = ci_95_u), linewidth = .5) +
#   geom_linerange(aes(xmin = ci_90_l, xmax = ci_90_u), linewidth = 1.25) +
#   geom_point(aes(x = estimate, fill = after_scale(clr_lighten(color))),
#              shape = 21, size = 2) +
#   facet_wrap(stat ~ ., scales = "free")

# (data_boots |> 
#     filter(demtype == "bot10-nes") |> 
#     pluck("bs_data"))[[1]] |> 
#   ggplot(aes(x = NCUR)) +
#   geom_density() +
#   geom_linerange(data = all_boots |> 
#                    select(demtype, stat, starts_with("ci_")) |>
#                    filter(stat == "NCUR", demtype == "bot10-nes"),
#                  aes(xmin = ci_95_l, xmax = ci_95_u, y = -1e-5), linewidth = .5, inherit.aes = FALSE) +
#   geom_linerange(data = all_boots |> select(demtype, stat, starts_with("ci_")) |> filter(stat == "NCUR", demtype == "bot10-nes"),
#                  aes(xmin = ci_90_l, xmax = ci_90_u, y = -1e-5), linewidth = 1.5, inherit.aes = FALSE) +
#   geom_point(aes(y = -2e-5), size = .3) +
#   geom_point(data = all_estimates |> filter(stat == "NCUR", demtype == "bot10-nes"),
#              aes(y = 0, x = estimate))


# CI based on 2.5 % and 97.5% percentiles of the bootstraped values
clr_ano <- clr_default[2]
all_estimates |> 
  left_join(data_boots |> select(demtype, data_quantiles) |> unnest(data_quantiles)) |> 
  ggplot(aes(y = demtype, color = demtype)) +
  geom_linerange(aes(xmin = ci_95_l, xmax = ci_95_u), linewidth = .5) +
  geom_linerange(aes(xmin = ci_90_l, xmax = ci_90_u), linewidth = 1.25) +
  # geom_linerange(aes(xmin = ci_95_l, xmax = ci_95_u), linewidth = .5) +
  # geom_linerange(aes(xmin = mean - sd, xmax = mean + sd), linewidth = 1.25) +
  geom_point(aes(x = estimate, fill = after_scale(clr_lighten(color))),
             shape = 21, size = 2) +
  guides(color = guide_legend(ncol = 3)) +
  facet_wrap(stat ~ ., scales = "free") +
  theme_linedraw(base_family = fnt_sel) +
  theme(strip.background = element_rect(fill = clr_ano, color = clr_ano),
        panel.background = element_rect(fill = "transparent", color = clr_ano),
        panel.grid = element_line(color = clr_ano, linetype = 3, linewidth = .5),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,0))

all_estimates |> 
  left_join(data_boots |> select(demtype, data_quantiles) |> unnest(data_quantiles)) |> 
  write_tsv("~/Dropbox/David/elephant_seals/data/mirang_fastsimcoal_estimates.tsv")

ggsave(here("results/img/demography/parameter_ci.pdf"),
       width = 16, height = 5, device = cairo_pdf)


# test <- data_boots |> select(demtype, data_quantiles) |> unnest(data_quantiles)
# test_q <- "bot06-lgm"
# test_s <- "NBOT"
# (data_boots |> 
#     filter(demtype == test_q) |> 
#     pluck("bs_data"))[[1]] |> 
#   ggplot(aes(x = .data[[test_s]])) +
#   geom_density() +
#   geom_linerange(data = test |> select(demtype, stat, starts_with("ci_")) |> filter(stat == test_s, demtype == test_q),
#                  aes(xmin = ci_95_l, xmax = ci_95_u, y = -1e-5), linewidth = .5, inherit.aes = FALSE) +
#   geom_linerange(data = test |> select(demtype, stat, starts_with("ci_")) |> filter(stat == test_s, demtype == test_q),
#                  aes(xmin = ci_90_l, xmax = ci_90_u, y = -1e-5), linewidth = 1.5, inherit.aes = FALSE) +
#   geom_jitter(aes(y = -2e-5), size = .5, color = clr_alpha("red"),width = .5, height = .01) +
#   geom_point(data = all_estimates |> filter(stat == test_s, demtype == test_q),
#              aes(y = 0, x = estimate), size = 3, shape = 1, color = "red")

# ================
#  comparison rad
# ================

estimates_rad <- tribble(
  ~stat,	~Estimate,	~lowerCI, ~upperCI,
  "NGM",	266.5,	226.5,	1721.5,
  "T1",	1222,	350,	1570,
  "NANC",	12855.5,	2827.5,	20274.5,
  "NBOT",	5.5,	5,	7.5,
  "NCUR",	2624,	2506.5,	2773 ) |> 
  mutate(demtype = "RAD")

all_estimates |> 
  left_join(data_boots |> select(demtype, data_quantiles) |> unnest(data_quantiles)) |> 
  mutate(across(c(estimate, mean:ci_95_u), ~ if_else(grepl("N", stat), .x / 2, .x) )) |> 
  ggplot(aes(y = demtype, color = demtype)) +
  geom_linerange(aes(xmin = ci_95_l, xmax = ci_95_u), linewidth = .5) +
  geom_linerange(aes(xmin = ci_90_l, xmax = ci_90_u), linewidth = 1.25) +
  geom_linerange(data = estimates_rad, aes(xmin = lowerCI, xmax = upperCI), linewidth = .5) +
  geom_point(data = estimates_rad, aes(x = Estimate, fill = after_scale(clr_lighten(color))),
             shape = 21, size = 2) +
  # geom_linerange(aes(xmin = ci_95_l, xmax = ci_95_u), linewidth = .5) +
  # geom_linerange(aes(xmin = mean - sd, xmax = mean + sd), linewidth = 1.25) +
  geom_point(aes(x = estimate, fill = after_scale(clr_lighten(color))),
             shape = 21, size = 2) +
  guides(color = guide_legend(ncol = 3)) +
  facet_wrap(stat ~ ., scales = "free") +
  theme_linedraw(base_family = fnt_sel) +
  theme(strip.background = element_rect(fill = clr_ano, color = clr_ano),
        panel.background = element_rect(fill = "transparent", color = clr_ano),
        panel.grid = element_line(color = clr_ano, linetype = 3, linewidth = .5),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,0))

ggsave(here("results/img/demography/parameter_ci_with_RAD.pdf"),
       width = 16, height = 5, device = cairo_pdf)

# ---------------
library(ggdist)
library(ggnewscale)
library(patchwork)
n_levels <- c("NCUR",  "NBOT", "NANC", "NLGM")

data_bot_expanded <- data_boots |>
  select(demtype, bs_data) |>
  unnest(bs_data) |> 
  select(demtype, NCUR, NANC, NBOT, NLGM = NGM, T1)

plt_bs_estimate <- \(param,
                     col = clr1,
                     col2 = clr2,
                     adj = .45,
                     n_fine = 301,
                     size_est = 3){
  data_est <- all_estimates |>
    filter(str_detect(demtype, "lgm")) |> 
    pivot_wider(values_from = estimate, names_from = stat) |> 
    rename(NLGM = "NGM") |> 
    mutate(across(starts_with("N"), \(x){x/2}),
           demtype = str_remove(demtype, "-lgm"))
  
  data_bot_expanded |>
    filter(str_detect(demtype, "lgm")) |> 
    mutate(across(starts_with("N"), \(x){x/2}),
           demtype = str_remove(demtype, "-lgm")) |> 
    ggplot(aes(y = demtype, x = .data[[param]])) +
    stat_slab(color = "transparent", fill = "transparent",n = 2) +
    geom_linerange(data = data_est,
                   linewidth = .3,
                   aes(color = demtype,
                       ymin = as.numeric(factor(demtype))-.45,
                       ymax = as.numeric(factor(demtype))+.45)) +
    stat_slab(data = ~ subset(., str_detect(demtype, "06")),
              aes(thickness = after_stat(pdf),
                  fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, q_quant)))),
              color = col, size = .25, scale = .75, side = "both",
              fill  = "white",
              adjust = adj,
              normalize = "groups",
              trim = FALSE, n = n_fine) +
    scale_colour_ramp_discrete(from = col, aesthetics = "fill_ramp", guide = "none") +
    new_scale(new_aes = "fill_ramp") +
    stat_slab(data = ~ subset(., !str_detect(demtype, "06")),
              aes(thickness = after_stat(pdf),
                  fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, q_quant)))),
              color = col2, size = .25, scale = .75, side = "both",
              fill  = "white",
              adjust = adj,
              normalize = "groups",
              trim = FALSE, n = n_fine) +
    scale_colour_ramp_discrete(from = col2, aesthetics = "fill_ramp", guide = "none") +
    stat_dots(side = "both", scale = 0.3, quantiles = 100,
              color = clr_alpha("black"), shape = 19) +
    geom_point(data = data_est,
               shape = 21,
               size = size_est,
               stroke = .4,
               aes(color = demtype), 
               fill = clr_alpha("white")) +
    scale_color_manual(values = c(`bot06` = clr_darken(col,.3), 
                                  `bot10` = clr_darken(col2,.3), 
                                  `null` = clr_darken(col2,.3)),
                       guide = "none") +
    labs(x = str_c(param, p_units[[param]])) +
    theme_bw(base_family = "Arial") +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          axis.title.y = element_blank())
}

p0 <- all_estimates |> 
  filter(demtype == "bot06-lgm") |> 
  left_join(data_boots |> select(demtype, data_quantiles) |> unnest(data_quantiles)) |>
  filter(grepl("^N", stat)) |> 
  mutate(across(c(estimate, mean:ci_95_u), \(x){x/2}),
         stat = ifelse(stat == "NGM", "NLGM", stat) |> factor(levels = n_levels)) |> 
  ggplot(aes(x = as.numeric(stat))) +
  geom_ribbon(aes(ymin = ci_95_l, max = ci_95_u), fill = clr_lighten(clr_accent,.45)) +
  geom_ribbon(aes(ymin = ci_90_l, max = ci_90_u), fill = clr_accent2) +
  geom_line(aes(y = estimate), color = "white", linetype = 3) +
  geom_point(aes(y = estimate),
             shape = 21,
             size = 3,
             stroke = .4,
             color = clr_accent, 
             fill = clr_alpha("white")) +
  scale_x_reverse(breaks = 4:1, labels = str_c(rev(n_levels), c("", " (T1)", "", ""))) +
  coord_cartesian(expand = 0,
                  xlim = c(4.02, .9),
                  ylim = c(-100, 14000)) +
  labs(x = "Timepoint", y = "Population Size (*N<sub>e</sub>*)") +
  theme_bw(base_family = "Arial") +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        axis.title.y = ggtext::element_markdown())

p_units <- c(NCUR = " (*N<sub>e</sub>*)",
             NBOT = " (*N<sub>e</sub>*)",
             NLGM = " (*N<sub>e</sub>*)",
             NANC = " (*N<sub>e</sub>*)",
             T1 = " (generations)")

p_list <- c(rev(n_levels), "T1") |> 
  map(plt_bs_estimate, adj = .75, col = clr_accent)

wrap_plots(c(p_list, list(p0)))  +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.margin = margin(r = unit(5,"pt")))) &
  theme(plot.tag = element_text(family = "Arial"),
        axis.title.x = ggtext::element_markdown(),
        panel.grid = element_blank())

ggsave("results/img/demography/param_selection.svg", width = 10, height = 5)
ggsave("results/img/demography/param_selection.pdf", width = 10, height = 5, device = cairo_pdf)

# main figure

plt_bs_estimate_single <- \(param,
                     col = clr1,
                     col2 = clr2,
                     adj = .45,
                     n_fine = 301,
                     size_est = 3){
  data_est <- all_estimates |>
    filter(str_detect(demtype, "lgm"),
           str_detect(demtype, "06")) |> 
    pivot_wider(values_from = estimate, names_from = stat) |> 
    rename(NLGM = "NGM") |> 
    mutate(across(starts_with("N"), \(x){x/2}),
           demtype = str_remove(demtype, "-lgm"))
  
  data_bot_expanded |>
    filter(str_detect(demtype, "lgm"),
           str_detect(demtype, "06")) |> 
    mutate(across(starts_with("N"), \(x){x/2}),
           demtype = str_remove(demtype, "-lgm")) |> 
    ggplot(aes(y = demtype, x = .data[[param]])) +
    stat_slab(color = "transparent", fill = "transparent",n = 2) +
    geom_linerange(data = data_est,
                   linewidth = .3,
                   aes(color = demtype,
                       ymin = as.numeric(factor(demtype))-.45,
                       ymax = as.numeric(factor(demtype))+.45)) +
    stat_slab(aes(thickness = after_stat(pdf),
                  fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, q_quant)))),
              color = col, size = .25, scale = .75, side = "both",
              fill  = "white",
              adjust = adj,
              normalize = "groups",
              trim = FALSE, n = n_fine) +
    scale_colour_ramp_discrete(from = col, aesthetics = "fill_ramp", guide = "none") +
    stat_dots(side = "both", scale = 0.3, quantiles = 100,
              color = clr_alpha("black"), shape = 19) +
    geom_point(data = data_est,
               shape = 21,
               size = size_est,
               stroke = .4,
               aes(color = demtype), 
               fill = clr_alpha("white")) +
    scale_color_manual(values = c(`bot06` = clr_darken(col,.3), 
                                  `bot10` = clr_darken(col2,.3), 
                                  `null` = clr_darken(col2,.3)),
                       guide = "none") +
    labs(x = "*N<sub>e</sub>*",
        y = param) +
    theme_bw(base_family = "Arial") +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          axis.title.y = ggtext::element_markdown(),
          # axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.subtitle = element_text(hjust = .5))
}

patch1 <- rev(n_levels) |> 
  map(plt_bs_estimate_single, adj = .75, col = clr_accent) |> 
  map2(.y = c(list(list(coord_flip(ylim = c(0.4,1.6)))),
              rep(list(list(coord_flip(ylim = c(0.4,1.6)),
                       theme(axis.title.y = element_blank()))),3)),
       .f = \(x,y){x + y}) |> 
  wrap_plots(nrow = 1)

p_t <- plt_bs_estimate_single("T1", adj = .75, col = clr_accent) +
  coord_cartesian(ylim = c(0.4,1.6)) +
  labs(y = "T1", x = "generations") +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(),
        plot.subtitle = element_text(hjust = .5))

((patch1 /
  p_t +
  plot_layout(heights = c(1,.25))) |
  p0) &
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.margin = margin(r = unit(5,"pt")))) &
  theme(plot.tag = element_text(family = "Arial"),
        panel.grid = element_blank())

ggsave("results/img/demography/param_selection_main.pdf", width = 9, height = 4.5, device = cairo_pdf)

