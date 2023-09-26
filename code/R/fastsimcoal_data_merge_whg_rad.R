library(tidyverse)
library(here)
library(glue)

## --- raw data for density ---
# rad functions
read_rad <- \(tag){
  read_tsv(here(glue("results/RAD/demography/Non_Param_boots_{tag}Model.txt"))) |> 
    mutate(demtype = rad_translate[tag])
}

# whg functions
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

get_boot_stats <- \(type){
  tibble(demtype = str_c("whg_", str_remove(type, "-lgm")),
         bs_data = list(0:99 |> 
                          map_dfr(read_bs, type = type) |> 
                          mutate(bs_dx = 0:99))
  )
}
# -----------
basepath <- "results/demography/fastsimcoal/mirang_on_mirang/"

rad_tags <- c("Bott10", "Bott6_Best", "Null")
rad_translate <- c("rad_bot10", "rad_bot06", "rad_null") |> 
  set_names(rad_tags)

data_boot_rad <- rad_tags |> 
  map_dfr(read_rad)

data_boot_rad_extra <- data_boot_rad |> 
  filter(demtype == "rad_bot06") |> 
  bind_rows(read_tsv(here("results/RAD/demography/NonParam_900ExtraBoots.txt"))) |> 
  mutate(demtype = "rad_bot06_1k")

all_dem_types <- c("bot06-lgm", "bot10-lgm", "null-lgm")

config_table <- tibble(type = all_dem_types) |> 
  mutate(stat = map(type, get_parameter_inventory)) |> 
  unnest(stat)

data_boot_whg <- all_dem_types |> 
  map_dfr(get_boot_stats) |> 
  unnest(bs_data) |> 
  select(-bs_dx) |> 
  mutate(across(starts_with("N"), \(x){x/2})) # haplotype ti individual

data_boot <- data_boot_rad |> 
  bind_rows(data_boot_rad_extra) |> 
  bind_rows(data_boot_whg) |>
  rename(NLGM = "NGM")

data_boot |> 
  write_tsv(here("results/demography/all_models_raw_boots.tsv"))

## --- confidence interval data ---
q_quant <- c(.95, 0.66)
q_borders <- c((1 - q_quant[1])/2, (1 - q_quant[2])/2)

data_ci <- data_boot |> 
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
  separate(q_in, into = c("ci_66_l", "ci_66_u"), sep = "_", convert = TRUE) |> 
  separate(q_out, into = c("ci_95_l", "ci_95_u"), sep = "_", convert = TRUE)

data_ci |> 
  write_tsv(here("results/demography/all_models_ci_boots.tsv"))

## --- point estimate data ----
# whg functions
get_estimates <- \(type){
  read_tsv(here(glue("{basepath}/{type}/bestrun/mirang_on_mirang_{type}.bestlhoods"))) |> 
    # pivot_longer(everything(), names_to = "stat", values_to = "estimate") |> 
    mutate(demtype = str_c("whg_", str_remove(type, "-lgm"))) |> 
    select(demtype, starts_with("N"), T1)
}
# -----------
rad_estimates <- read_tsv(here("results/RAD/demography/FSC_parameterEstiamtes.csv"),
                                  skip = 1) |> 
  select(-`...2`)  |> 
  mutate(n = c(1,2,1)) |> 
  uncount(n) |> 
  mutate(
    Model = if_else(duplicated(Model), "Bott_6_1k", Model),
    demtype = case_when(
    Model == "Null" ~ "rad_null",
    Model == "Bott_6" ~ "rad_bot06",
    Model == "Bott_6_1k" ~ "rad_bot06_1k",
    Model == "Bott_10" ~ "rad_bot10",
    .default = NA
  ))  |> 
  select(demtype, starts_with("N"), T1)

whg_estimates <- all_dem_types |> 
  map_dfr(get_estimates) |> 
  mutate(across(starts_with("N"), \(x){x/2})) # haplotype ti individual

data_estimates <- rad_estimates |> 
  bind_rows(whg_estimates) |>
  rename(NLGM = "NGM")

data_estimates |> 
  write_tsv(here("results/demography/all_models_estimates.tsv"))
