library(tidyverse)
library(prismatic)
library(here)
library(glue)
library(boot)
source(here("code/R/project_defaults.R"))
clr1 <- "red"
clr2 <- "gray20"
basepath <- "results/demography/fastsimcoal/mirang_on_mirang/"

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
                                conf = c(0.90, 0.95),
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
                                   q90 = str_c(quantile(x, c(.05, .95), na.rm = TRUE), collapse = "_"),
                                   q95 = str_c(quantile(x, c(.025, .975), na.rm = TRUE), collapse = "_")))
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
         bs_data = list(1:50 |> 
                          map_dfr(read_bs, type = type) |> 
                          mutate(bs_dx = 1:50))
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
clr_ano <- "gray70"
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
