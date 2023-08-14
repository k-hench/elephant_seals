library(tidyverse)
library(prismatic)
library(here)
library(ggtext)
library(glue)
library(patchwork)
library(ggstance)
source(here("code/R/project_defaults.R"))

c_select <- rcartocolor::carto_pal(12, "Prism") |>  color()
clr_sfs <- c_select[c(1, 4, 2, 8, 9)] |>  set_names(nm = c("whg", "RAD", "RAD_binned", "leo", "whg_binned"))

basepath <- "results/demography/mirang_on_mirang/"

import_sfs <- \(type, demtype, prefix){
  if(type == "obs"){
    n_skip <- 1
    fl <- here(str_c("results/demography/sfs/mirang_on_mirang/fastsimcoal2/mirang_MAFpop0.obs"))
  } else if(type == "leo") {
    n_skip <- 1
    fl <- "~/Downloads/sfs/mirleo_MAFpop0.obs"
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
sfs_leo <- "~/Downloads/sfs/mirleo_MAFpop0.obs"
sfs_whg <- here(str_c("results/demography/sfs/mirang_on_mirang/fastsimcoal2/mirang_MAFpop0.obs"))

data_rad <- tibble(snp_count =  str_split(read_lines(sfs_rad), pattern = " ")[[1]] |> as.numeric()) |> 
  mutate(allele_count = seq_along(snp_count) -1,
         sfs_type = "RAD",
         allele_freq = allele_count / max(allele_count),
         bin_width = 1 / max(allele_count)) |>
    filter(allele_count != 0) |> 
  mutate(snp_freq = snp_count / sum(snp_count)) |> 
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

data_rad_binned_20 <- tibble(snp_count =  str_split(read_lines(sfs_rad),
                                                 pattern = " ")[[1]] |> as.numeric()) |> 
  mutate(allele_count = seq_along(snp_count) -1,
         allele_freq = allele_count / max(allele_count),
         allele_count_bin = cut(allele_freq, breaks = (0:20)/20)) |>
  filter(allele_count != 0) |> 
  group_by(allele_count_bin) |> 
  summarise(allele_count = max(allele_count),
            allele_freq = max(allele_freq),
            snp_count = sum(snp_count)) |> 
  ungroup() |>
  select(-allele_count_bin) |> 
  mutate(sfs_type = "RAD",
         snp_freq = snp_count / sum(snp_count),
         bin_width = 1 / 21) |> 
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
         bin_width = 1 / 40)

data_whg_binned_20 <- read_lines(sfs_whg, skip = 1, n_max = 2) |>
  str_replace_all("\t", " ") |>
  str_remove_all("d0_") |>
  str_split(" ") |>
  set_names(nm = c("allele_count", "snp_count")) |>
  as_tibble() |>
  mutate(across(everything(), as.numeric)) |> 
  filter(!is.na(allele_count)) |> 
  mutate(allele_count = seq_along(snp_count) -1,
         allele_freq = allele_count / max(allele_count),
         allele_count_bin = cut(allele_freq, breaks = (0:20)/20)) |>
  filter(allele_count != 0) |> 
  group_by(allele_count_bin) |> 
  summarise(allele_count = max(allele_count),
            allele_freq = max(allele_freq),
            snp_count = sum(snp_count)) |> 
  ungroup() |>
  select(-allele_count_bin) |> 
  mutate(sfs_type = "whg",
         snp_freq = snp_count / sum(snp_count),
         bin_width = 1 / 21) |> 
  select(allele_count, snp_count, sfs_type, snp_freq, allele_freq, bin_width)

data_leo <- read_lines(sfs_leo, skip = 1, n_max = 2) |>
  str_replace_all("\t", " ") |>
  str_remove_all("d0_") |>
  str_split(" ") |>
  set_names(nm = c("allele_count", "snp_count")) |>
  as_tibble() |>
  mutate(across(everything(), as.numeric)) |> 
  filter(!is.na(allele_count)) |> 
  mutate(sfs_type = "leo",
         snp_freq = snp_count / sum(snp_count),
         allele_freq = allele_count / max(allele_count),
         bin_width = 1 / 20)

ggplot(mapping = aes(x = allele_freq - bin_width)) +
  geom_step(data = data_whg, aes(y = snp_freq, color = "whg")) +
  geom_step(data = data_leo, aes(y = snp_freq, color = "leo")) +
  geom_step(data = data_rad, aes(y = snp_freq * 1.5e3, color = "RAD")) +
  geom_step(data = data_rad_binned, aes(y = snp_freq * 4e2, color = "RAD_binned")) +
  coord_cartesian(ylim = c(0,.3)) +
  scale_color_manual("SFS type",
                     values = clr_sfs) +
  theme_minimal(base_family = fnt_sel) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1))

ggplot(mapping = aes(x = allele_freq - bin_width)) +
  geom_step(data = data_whg, aes(y = snp_count, color = "whg")) +
  geom_step(data = data_leo, aes(y = snp_count, color = "leo")) +
  geom_step(data = data_rad, aes(y = snp_count * 1.5e3, color = "RAD")) +
  geom_step(data = data_rad_binned, aes(y = snp_count * 4e2, color = "RAD_binned")) +
  # coord_cartesian(ylim = c(0,1e5)) +
  facet_grid(sfs_type ~ ., scales = "free_y") +
  scale_color_manual("SFS type",
                     values = clr_sfs) +
  theme_minimal(base_family = fnt_sel) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1))

ggsave(filename = here("~/Dropbox/David/elephant_seals/img/demography/compare_rad_whg_SFS_binned_only.pdf"),
       width = 6, height = 3.5, device = cairo_pdf)

ggplot(mapping = aes(x = allele_freq - bin_width)) +
  geom_step(data = data_whg_binned_20, aes(y = snp_freq, color = "whg_binned")) +
  geom_step(data = data_whg, aes(y = snp_freq, color = "whg")) +
  geom_step(data = data_leo, aes(y = snp_freq, color = "leo")) +
  geom_step(data = data_rad_binned_20, aes(y = snp_freq, color = "RAD_binned")) +
  geom_step(data = data_rad, aes(y = snp_freq, color = "RAD")) +
  coord_cartesian(ylim = c(0,1)) +
  facet_grid(sfs_type ~ ., scales = "free_y") +
  scale_color_manual("SFS type",
                     values = clr_sfs) +
  theme_minimal(base_family = fnt_sel) +
  theme(legend.position = "bottom")

ggsave(filename = "~/Dropbox/David/elephant_seals/img/demography/compare_rad_whg_leo_SFS_binned_only.pdf",
       width = 6, height = 5, device = cairo_pdf)
