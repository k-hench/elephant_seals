library(tidyverse)
library(prismatic)
library(here)
library(glue)
source(here("code/R/project_defaults_shared.R"))

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
                                   q_in = str_c(quantile(x, c(1 - q_borders[2], q_borders[2]), na.rm = TRUE), collapse = "_"),
                                   q_out = str_c(quantile(x, c(1 - q_borders[1], q_borders[1]), na.rm = TRUE), collapse = "_")))
                     })) |> 
    ungroup() |> 
    pivot_longer(everything(), names_to = "stat") |> 
    unnest(value) |> 
    separate(q_in, into = c("ci_in_l", "ci_in_u"), sep = "_", convert = TRUE) |> 
    separate(q_out, into = c("ci_out_l", "ci_out_u"), sep = "_", convert = TRUE)
}

all_dem_types <- c("bot06-lgm", "bot06-nes",
                   "bot10-lgm", "bot10-nes",
                   "null-lgm", "null-nes")

config_table <- tibble(type = all_dem_types) |> 
  mutate(stat = map(type, get_parameter_inventory)) |> 
  unnest(stat)

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

data_boots |> 
  select(demtype, data_quantiles) |> 
  unnest(data_quantiles) |> 
  write_tsv(here("results/demography/fastsimcoal/mirang_on_mirang/bootstrap_ci.tsv.gz"))

data_boots |>
  select(demtype, bs_data) |>
  unnest(bs_data) |> 
  select(demtype, NCUR, NANC, NBOT, NLGM = NGM, T1) |> 
  write_tsv(here("results/demography/fastsimcoal/mirang_on_mirang/bootstrap_ci_expanded.tsv.gz"))
