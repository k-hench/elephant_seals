library(tidyverse)
library(prismatic)
library(here)
library(ggtext)
library(glue)
library(patchwork)
library(ggstance)
source(here("code/R/project_defaults.R"))
clr1 <- "red"
clr2 <- "gray20"
basepath <- "results/demography/mirang_on_mirang/"

import_sfs <- \(type, demtype, prefix){
  if(type == "obs"){
    n_skip <- 1
    fl <- here(str_c("results/demography/sfs/mirang_on_mirang/fastsimcoal2/mirang_MAFpop0.obs"))
  } else if(type == "exp") {
    n_skip <- 0
    fl <- here(str_c(basepath, demtype,"/bestrun/",prefix,"_MAFpop0.txt"))
  }
  read_lines(fl, skip = n_skip, n_max = 2) |>
    str_replace_all("\t", " ") |>
    str_remove_all("d0_") |>
    str_split(" ") |>
    set_names(nm = c("allele_count", "snp_count")) |>
    as_tibble() |>
    mutate(across(everything(), as.numeric)) |> 
    filter(!is.na(allele_count)) |> 
    mutate(sfs_type = type,
           snp_freq = snp_count / sum(snp_count))
}

sfs_rad <- "~/Dropbox/David/elephant_seals/img/demography/RADdata_withXchr_SFS.sfs"
sfs_whg <- here(str_c("results/demography/sfs/mirang_on_mirang/fastsimcoal2/mirang_MAFpop0.obs"))

data_rad <- tibble(snp_count =  str_split(read_lines(sfs_rad), pattern = " ")[[1]] |> as.numeric()) |> 
  mutate(allele_count = seq_along(snp_count) -1,
         sfs_type = "RAD",
         snp_freq = snp_count / sum(snp_count),
         allele_freq = allele_count / max(allele_count),
         bin_width = 1 / max(allele_count)) |>
    filter(allele_count != 0) |> 
  select(allele_count, snp_count, sfs_type, snp_freq, allele_freq, bin_width)


data_rad_binned <- tibble(snp_count =  str_split(read_lines(sfs_rad),
                                                 pattern = " ")[[1]] |> as.numeric()) |> 
  mutate(allele_count = seq_along(snp_count) -1,
         allele_freq = allele_count / max(allele_count),
         allele_count_bin = cut(allele_freq, breaks = (0:40)/40)) |>
  filter(allele_count != 0) |> 
  group_by(allele_count_bin) |> 
  summarise(allele_count = max(allele_count),
            allele_freq = max(allele_freq),
            snp_count = sum(snp_count)) |> 
  ungroup() |>
  select(-allele_count_bin) |> 
  mutate(sfs_type = "RAD",
         snp_freq = snp_count / sum(snp_count),
         bin_width = 1 / 41) |> 
  select(allele_count, snp_count, sfs_type, snp_freq, allele_freq, bin_width)

data_whg <- read_lines(sfs_whg, skip = 1, n_max = 2) |>
  str_replace_all("\t", " ") |>
  str_remove_all("d0_") |>
  str_split(" ") |>
  set_names(nm = c("allele_count", "snp_count")) |>
  as_tibble() |>
  mutate(across(everything(), as.numeric)) |> 
  filter(!is.na(allele_count)) |> 
  mutate(sfs_type = "whg",
         snp_freq = snp_count / sum(snp_count),
         allele_freq = allele_count / max(allele_count),
         bin_width = 1 / max(allele_count))

ggplot(mapping = aes(x = allele_freq - .5 * bin_width)) +
  geom_step(data = data_whg, aes(y = snp_count, color = "whg")) +
  geom_step(data = data_rad, aes(y = snp_count * 1.5e3, color = "RAD")) +
  geom_step(data = data_rad_binned, aes(y = snp_count * 4e2, color = "RAD_binned")) +
  coord_cartesian(ylim = c(0,1e5)) +
  scale_color_manual("SFS type",
                     values = c(RAD = clrs[["mirleo"]], RAD_binned = clr_darken(clrs[["mirleo"]],.3), whg = "black")) +
  theme_minimal(base_family = fnt_sel) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1))


ggsave(filename = here("~/Dropbox/David/elephant_seals/img/demography/compare_rad_whg_SFS_binned_only.pdf"),
       width = 6, height = 3.5, device = cairo_pdf)

