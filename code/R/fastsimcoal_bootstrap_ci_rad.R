library(tidyverse)
library(prismatic)
library(here)
library(glue)
source(here("code/R/project_defaults_shared.R"))

basepath <- "results/demography/fastsimcoal/mirang_on_mirang/"

get_estimates <- \(type){
  read_tsv(here(glue("{basepath}/{type}/bestrun/mirang_on_mirang_{type}.bestlhoods"))) |> 
    pivot_longer(everything(), names_to = "stat", values_to = "estimate") |> 
    mutate(demtype = type)
}

all_dem_types <- c("bot06-lgm", "bot06-nes",
                   "bot10-lgm", "bot10-nes",
                   "null-lgm", "null-nes")

all_estimates <- all_dem_types |> 
  map_dfr(get_estimates)

data_boots <- read_tsv(here("results/demography/fastsimcoal/mirang_on_mirang/bootstrap_ci.tsv.gz"))
clr_ano <- "lightgray"
# CI based on 2.5 % and 97.5% percentiles of the bootstraped values
all_estimates |> 
  left_join(data_boots) |> 
  ggplot(aes(y = demtype, color = demtype)) +
  geom_linerange(aes(xmin = ci_out_l, xmax = ci_out_u), linewidth = .5) +
  geom_linerange(aes(xmin = ci_in_l, xmax = ci_in_u), linewidth = 1.25) +
  geom_point(aes(x = estimate, fill = after_scale(clr_lighten(color))),
             shape = 21, size = 2) +
  guides(color = guide_legend(ncol = 3)) +
  facet_wrap(stat ~ ., scales = "free") +
  theme_ms()+
  # theme_linedraw(base_family = fnt_sel) +
  theme(#strip.background = element_rect(fill = clr_ano, color = clr_ano),
        panel.background = element_rect(fill = "transparent", color = clr_ano),
        #panel.grid = element_line(color = clr_ano, linetype = 3, linewidth = .5),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,0))

# ================
#  comparison rad
# ================
q_quant <- c(.95, 0.66)
q_borders <- c((1 - q_quant[1])/2, (1 - q_quant[2])/2)

rad_tags <- c("Bott10", "Bott6_Best", "Null")
rad_translate <- str_c("RAD_",
                       c("bot10", "bot06", "null")) |> 
  str_c("-lgm") |> 
  set_names(rad_tags)

read_rad <- \(tag){
  read_tsv(here(glue("results/RAD/demography/Non_Param_boots_{tag}Model.txt"))) |> 
    mutate(demtype = rad_translate[tag])
}

data_boot_rad <- rad_tags |> 
  map_dfr(read_rad)

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

estimates_rad <- tribble(
  ~stat,	~Estimate,	~lowerCI, ~upperCI,
  "NGM",	266.5,	226.5,	1721.5,
  "T1",	1222,	350,	1570,
  "NANC",	12855.5,	2827.5,	20274.5,
  "NBOT",	5.5,	5,	7.5,
  "NCUR",	2624,	2506.5,	2773 ) |> 
  mutate(demtype = "RAD_bot06-lgm")

all_estimates |> 
  left_join(data_boots) |> 
  mutate(across(c(estimate, mean:ci_out_u), ~ if_else(grepl("N", stat), .x / 2, .x) )) |> 
  ggplot(aes(y = demtype, color = demtype)) +
  geom_linerange(aes(xmin = ci_out_l, xmax = ci_out_u), linewidth = .5) +
  geom_linerange(aes(xmin = ci_in_l, xmax = ci_in_u), linewidth = 1.25) +
  # geom_linerange(data = estimates_rad, aes(xmin = lowerCI, xmax = upperCI), linewidth = .5) +
  geom_linerange(data = data_boot_ci,
                 aes(xmin = ci_out_l, xmax = ci_out_u), linewidth = .5) +
  geom_linerange(data = data_boot_ci,
                 aes(xmin = ci_in_l, xmax = ci_in_u), linewidth = 1.25) +
  geom_point(data = estimates_rad, aes(x = Estimate, fill = after_scale(clr_lighten(color))),
             shape = 21, size = 2) +
  geom_point(aes(x = estimate, fill = after_scale(clr_lighten(color))),
             shape = 21, size = 2) +
  guides(color = guide_legend(ncol = 3)) +
  facet_wrap(stat ~ ., scales = "free")  +
  theme_ms()+
  # theme_linedraw(base_family = fnt_sel) +
  theme(#strip.background = element_rect(fill = clr_ano, color = clr_ano),
    panel.background = element_rect(fill = "transparent", color = clr_ano),
    #panel.grid = element_line(color = clr_ano, linetype = 3, linewidth = .5),
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

data_bot_expanded <- read_tsv(here("results/demography/fastsimcoal/mirang_on_mirang/bootstrap_ci_expanded.tsv.gz"))

estimates_rad_all_mod <- read_tsv(here("results/RAD/demography/FSC_parameterEstiamtes.csv"),
                                  skip = 1) |> 
  select(-`...2`) |>
  rename(NLGM = "NGM") |> 
  mutate(demtype = case_when(
    Model == "Null" ~ "rad-null",
    Model == "Bott_6" ~ "rad-bot06",
    Model == "Bott_10" ~ "rad-bot10",
    .default = NA
  )) |> 
  select(-Model) |> 
  select(demtype, everything())

plt_bs_estimate <- \(param,
                     col = clr_default[1],
                     col2 = clr_default[2],
                     adj = .45,
                     n_fine = 301,
                     size_est = 3,
                     target_dem = "rad-bot06"){
  data_est_whg <- all_estimates |>
    filter(str_detect(demtype, "lgm")) |> 
    pivot_wider(values_from = estimate, names_from = stat) |> 
    rename(NLGM = "NGM") |> 
    mutate(across(starts_with("N"), \(x){x/2}),
           demtype = str_c("whg-", demtype) |> str_remove("-lgm"))
  
  # data_est_rad <- estimates_rad |> 
  #   select(demtype, stat, Estimate) |> 
  #   pivot_wider(names_from = stat, values_from = Estimate) |> 
  #   mutate(demtype = str_replace(demtype, "RAD_", "rad-") |> 
  #            str_remove("-lgm")) |> 
  #   rename(NLGM = "NGM")
  
  data_est <- bind_rows(data_est_whg,
                        estimates_rad_all_mod)
  
  data_dist_whg <- data_bot_expanded |>
    filter(str_detect(demtype, "lgm")) |> 
    mutate(across(starts_with("N"), \(x){x/2}),
           demtype = str_c("whg-", demtype) |> str_remove("-lgm"))
  
  data_dist_rad <- data_boot_rad |> 
    select(demtype, starts_with("N"), T1) |> 
    mutate(demtype = str_replace(demtype, "RAD_", "rad-") |> 
             str_remove("-lgm"))|> 
    rename(NLGM = "NGM")
  
  data_dist <- bind_rows(data_dist_whg,
            data_dist_rad)
  data_dist |> 
    ggplot(aes(y = demtype, x = .data[[param]])) +
    stat_slab(color = "transparent", fill = "transparent", n = 2) +
    geom_linerange(data = data_est,
                   linewidth = .3,
                   aes(color = demtype,
                       ymin = as.numeric(factor(demtype, levels = sort(unique(data_dist$demtype)))) -.45,
                       ymax = as.numeric(factor(demtype, levels = sort(unique(data_dist$demtype)))) +.45)) +
    stat_slab(aes(thickness = after_stat(pdf),
                  fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, q_quant))),
                  fill  = demtype == target_dem,
                  color = demtype == target_dem),
              linewidth = .45, scale = .75, side = "both",
              adjust = adj,
              normalize = "groups",
              trim = FALSE, n = n_fine) +
    scale_colour_ramp_discrete(from = "white", aesthetics = "fill_ramp", guide = "none",range = c(1, 0)) +
    scale_fill_manual(values = c(`TRUE` = clr_default[[1]],
                                `FALSE` = clr_default[[2]]),
                    guide = "none") +
    scale_color_manual(values = c(`TRUE` = clr_darken(clr_default[[1]]),
                                  `FALSE` = clr_darken(clr_default[[2]])),
                       guide = "none") +
    stat_dots(side = "both", scale = 0.3, quantiles = 100,
              color = clr_alpha("black"), shape = 19) +
    geom_point(data = data_est,
               shape = 21,
               size = size_est,
               stroke = .4,
               aes(color = demtype), 
               fill = clr_alpha("white")) +
    labs(x = str_c(param, p_units[[param]])) +
    theme_ms() +
    theme(panel.border = element_blank(),
          axis.line = element_line(),
          axis.title.y = element_blank())
}

p0 <- estimates_rad |> 
  filter(demtype == "RAD_bot06-lgm") |> 
  select(stat, estimate = Estimate, demtype)|> 
  left_join(data_boot_ci) |>
  filter(grepl("^N", stat)) |> 
  mutate(stat = ifelse(stat == "NGM", "NLGM", stat) |> factor(levels = n_levels)) |> 
  ggplot(aes(x = as.numeric(stat))) +
  geom_ribbon(aes(ymin = ci_out_l, max = ci_out_u), fill = clr_lighten(clr_default[1],.45)) +
  geom_ribbon(aes(ymin = ci_in_l, max = ci_in_u), fill = clr_default[1]) +
  geom_line(aes(y = estimate), color = "white", linetype = 3) +
  geom_point(aes(y = estimate),
             shape = 21,
             size = 3,
             stroke = .4,
             color = clr_default[1], 
             fill = clr_alpha("white")) +
  scale_x_reverse(breaks = 4:1, labels = str_c(rev(n_levels), c("", " (T1)", "", ""))) +
  coord_cartesian(expand = 0,
                  xlim = c(4.02, .9),
                  ylim = c(-100, 20500)) +
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

ggsave("results/img/demography/param_selection_rad.svg", width = 10, height = 5)
ggsave("results/img/demography/param_selection_rad.pdf", width = 10, height = 5, device = cairo_pdf)

# main figure
plt_bs_estimate_single <- \(param,
                     col = clr1,
                     col2 = clr2,
                     adj = .45,
                     n_fine = 301,
                     size_est = 3){
  data_est <- estimates_rad |>
    filter(str_detect(demtype, "lgm"),
           str_detect(demtype, "06")) |> 
    select(demtype, stat, Estimate) |> 
    pivot_wider(names_from = stat, values_from = Estimate) |> 
    mutate(demtype = str_replace(demtype, "RAD_", "rad-") |> 
             str_remove("-lgm")) |> 
    rename(NLGM = "NGM")
    
  data_boot_rad  |>
    filter(str_detect(demtype, "lgm"),
           str_detect(demtype, "06")) |> 
    select(demtype, starts_with("N"), T1) |> 
    mutate(demtype = str_replace(demtype, "RAD_", "rad-") |> 
             str_remove("-lgm"))|> 
    rename(NLGM = "NGM") |> 
    ggplot(aes(y = demtype, x = .data[[param]])) +
    stat_slab(color = "transparent", fill = "transparent",n = 2) +
    geom_linerange(data = data_est,
                   linewidth = .3,
                   color = clr_default[[1]],
                   aes(ymin = as.numeric(factor(demtype))-.45,
                       ymax = as.numeric(factor(demtype))+.45)) +
    stat_slab(aes(thickness = after_stat(pdf),
                  fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, q_quant)))),
              color = clr_default[[1]],
              fill  = clr_default[[1]],
              linewidth = .45, scale = .75, side = "both",
              adjust = adj,
              normalize = "groups",
              trim = FALSE, n = n_fine) +
    scale_colour_ramp_discrete(from = "white", aesthetics = "fill_ramp", guide = "none", range = c(1,0)) +
    stat_dots(side = "both", scale = 0.3, quantiles = 100,
              color = clr_alpha("black"), shape = 19) +
    geom_point(data = data_est,
               shape = 21,
               size = size_est,
               stroke = .4,
               color = clr_darken(clr_default[[1]]), 
               fill = clr_alpha("white")) +
    labs(x = "*N<sub>e</sub>*",
        y = param) +
    theme_ms() +
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
  labs(y = "T1", x = "Generations") +
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
  p0 + theme(axis.line = element_line(linewidth = .3))) &
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.margin = margin(r = unit(5,"pt")))) &
  theme(plot.tag = element_text(family = "Arial"),
        panel.grid = element_blank())

ggsave("results/img/demography/param_selection_main_rad.pdf", width = 9, height = 4.5, device = cairo_pdf)

