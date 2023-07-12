library(tidyverse)
library(prismatic)
library(here)
library(ggtext)
source(here("code/R/project_defaults.R"))

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

get_both_sfs <- \(demtype){
  prefix <- str_c("mirang_on_mirang_", demtype)
  
  c("obs", "exp") |> 
    map_dfr(import_sfs, demtype = demtype, prefix = prefix) |> 
    mutate(demtype = demtype)
}

all_dem_types <- c("bot06-lgm", "bot06-nes",
                   "bot10-lgm", "bot10-nes",
                   "null-lgm", "null-nes")

data <- all_dem_types |> 
  map_dfr(get_both_sfs)

split_str_nth <- \(str, n){
  strsplit(str, sprintf("(?:[^,]*(?:,[^,]*){%s})\\K,", n-1), perl=TRUE)[[1]]
}

get_model_stats <- \(type){
  all_stats <- read_tsv(here(str_c("results/demography/mirang_on_mirang/", type, "/bestrun/mirang_on_mirang_", type, ".AIC")))[,2:1] |> 
    bind_cols(read_tsv(here(str_c("results/demography/mirang_on_mirang/", type, "/bestrun/mirang_on_mirang_", type, ".bestlhoods"))))
  all_labs <- str_c(str_c("<span style='color:#666666'>", names(all_stats), ":</span> "),
        round(all_stats[1,], digits = 2), collapse = ",") |> 
    split_str_nth(3) |> 
    str_c("<br>", collapse = "") |> 
    str_replace_all(",", ", ")
    
  tibble(demtype = type,
         label = all_labs)
}

parameter_estimates <- all_dem_types |> 
  map_dfr(get_model_stats)

p <- data |> 
  ggplot(aes(x = allele_count, y = snp_freq)) +
  geom_hline(yintercept = 0, linewidth = .1, color = "black") +
  geom_step(aes(color = sfs_type)) +
  geom_bar(data = data |>
             select(-snp_count) |>
             pivot_wider(values_from = snp_freq, names_from = sfs_type) |>
             mutate(diff = obs - exp),
           aes(x = allele_count + .5, y = diff, fill = "\U0394 (obs - exp)"),
           stat = "identity",
           color = "transparent") +
  geom_richtext(data = parameter_estimates,
                aes(x = 20, y = -Inf, label = label),
                vjust = .1, size = 4, family = fnt_sel,#"ubuntu mono",
                fill = NA, label.color = NA, # remove background and outline
                label.padding = grid::unit(rep(0, 4), "pt")) +
  scale_color_manual(NULL, values = c(exp = clr_alpha("red", .75),
                                      obs = clr_alpha("gray20", .75)),
                     labels = c(obs = "observed SFS", exp = "expected SFS")) +
  scale_fill_manual(NULL, values = clr_alpha("gray75")) +
  facet_wrap(demtype ~ . ) +
  labs(title = "Comparison of observed and expected folded 1D SFS for different demographic models",
       caption = "parameter estimates rounded to 2 digits",
       x = "Allele count", y = "SNP frequency") +
  theme_minimal(base_family = fnt_sel) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = .5),
        panel.grid = element_blank(),
        panel.background = element_rect(linewidth = .3, color = "gray50"))

ggsave(filename = here("results/img/demography/compare_exp_obs_SFS.pdf"),
       plot = p, width = 16, height = 9, device = cairo_pdf)
