library(tidyverse)
library(here)
library(glue)
library(readxl)

upper_3 <- \(x){str_replace(str_replace(x, "^whg", "wgs"), "^[a-z]{3}", toupper)}
get_lhoods <- \(type){
  read_tsv(here(glue("results/demography/fastsimcoal/mirang_on_mirang/{type}/bestrun/mirang_on_mirang_{type}.bestlhoods")),
           col_types = "dddddddddd") |> 
    mutate(demtype = str_c("WGS_", str_remove(type, "-lgm")))
}

ks <- c(bot06 = 5, bot10 = 5, null = 3)
params <- c(NLGM = "NeLGM", NANC = "NePREBOT", NBOT = "NeBOT", NCUR = "NePOSTBOT", T1 = "Tse")
all_dem_types <- c("bot06-lgm", "bot10-lgm", "null-lgm")

data_aic <- all_dem_types |>  
  map_dfr(\(type){ 
    read_tsv(here(str_c("results/demography/fastsimcoal/mirang_on_mirang/", type, "/bestrun/mirang_on_mirang_", type, ".AIC"))) |> 
      mutate(demtype = str_c("WGS_", str_remove(type, "-lgm")))}) |> 
  select(demtype , AIC)

data_lhood <- all_dem_types |>  
  map_dfr(get_lhoods) |> 
  select(demtype, MaxEstLhood, MaxObsLhood)

data_rad <- read_xlsx(here("results/RAD/demography/Supplementary_Table_David_1.xlsx")) |> 
  select(Model, MaxEstLhood = Lhood, AIC) |> 
  mutate(demtype = str_c("RAD_", c(Null = "null", Bott_6 = "bot06", Bott_10 = "bot10")[Model]),
         MaxObsLhood = -16686.145) |> 
  select(-Model)

data_estimates <- read_tsv(here("results/demography/all_models_estimates.tsv")) |> 
  mutate(demtype = upper_3(demtype)) |> 
  filter(demtype != upper_3("rad_bot06_1k")) |> 
  pivot_longer(-demtype, values_to = "estimate", values_drop_na = TRUE, names_to = "stat")

data_ci <- read_tsv(here("results/demography/all_models_ci_boots.tsv")) |> 
  mutate(demtype = upper_3(demtype)) |> 
  filter(demtype != upper_3("rad_bot06_1k")) |> 
  left_join(data_estimates, by = c("demtype", "stat")) |> 
  mutate(ci = str_c(sprintf('%.0f', estimate), " (",
                    sprintf('%.0f',median), ", [",
                    sprintf('%.0f',ci_95_u), ", ", 
                    sprintf('%.0f',ci_95_l), "])"),
         name = params[stat]) |> 
  select(demtype, name, ci) |>
  filter(!is.na(name)) |> 
  pivot_wider(names_from = name, values_from = ci) |> 
  mutate(across(everything(), \(x){str_replace(x, "NA \\(NA, \\[NA, NA\\]\\)", "-")}))

data_aic |> 
  left_join(data_lhood) |>
  bind_rows(data_rad) |> 
  left_join(data_ci) |> 
  separate(demtype, into = c("method", "model"), sep = "_") |> 
  mutate(k = ks[model],
         across(c(AIC, MaxEstLhood, MaxObsLhood), \(x){sprintf('%.0f', x)})) |> 
  arrange(method, model) |> 
  select(Method = method,
         Model = model,
         k,
         AIC,
         est_Lhood = MaxEstLhood,
         obs_Lhood = MaxObsLhood,
         NeLGM,
         NePREBOT,
         NeBOT,
         NePOSTBOT,
         Tse) |> 
  write_tsv(here("results/tab/suppl_model_summary.tsv"))
