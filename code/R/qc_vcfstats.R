library(tidyverse)
library(here)
library(ggstance)
library(prismatic)
source("code/R/project_defaults.R")

files <- c(here("results/genotyping/raw/mirang_raw_snps_stat.txt"),
           dir(here("results/genotyping/filtered/"), pattern = "txt",full.names = TRUE))

get_vcfstats <- \(file){
  tibble(raw = read_lines(file) |> 
    str_remove("kept ") |> 
    str_replace(" out of [a-z ]*", ",") |> 
    str_replace(" ", ",")) |> 
    separate(raw, into = c("n","n_total", "stat"), convert = TRUE) |> 
    select(-n_total) |> 
    # pivot_wider(names_from = stat,values_from = n) |> 
    mutate(file = str_remove(file, ".*/") |> str_remove("_stat.txt"),
           n = if_else(stat == "Sites", n*1e-6, n),
           stat = if_else(stat == "Sites", "SNPs (n X 10<sup>6</sup>)", "n Individuals"))
  }

p <- files |> 
  map_dfr(get_vcfstats) |> 
  mutate(file = factor(file, levels = c("mirang_raw_snps", "mirang_filtered_all",
                                        "mirang_bi-allelic", "mirang_mac1",
                                        "mirang_filtered_mirang", "mirang_filtered_mirleo") |> 
                         rev())) |> 
  ggplot(aes(x = n, y = file)) +
  geom_barh(stat = "identity",
            color = "gray60", fill = clr_alpha("gray85")) +
  facet_wrap(stat ~. , scales = "free", nrow = 2) +
  theme_minimal(base_family = fnt_sel) +
  theme(strip.text = ggtext::element_markdown())
  
ggsave(filename = here("results/img/qc/vcf_stats.pdf"),
       plot = p,
       width = 5, height = 5, device = cairo_pdf)
