library(tidyverse)
library(here)

data_het <- read_tsv(here("results/het/win_het_ind_all_w1Mb_s250kb.tsv.gz")) |> 
  group_by(spec, gpos) |> 
  filter(row_number() == 1) |> 
  ungroup() |> 
  arrange(gpos)

data_het |> 
  group_by(spec) |> 
  summarise(n_snps_mean = mean(n_snps),
            n_snps_median = median(n_snps),
            n_snps_sd = sd(n_snps),
            n_snps_min = mean(n_snps),
            n_snps_max = mean(n_snps),
            n_snps_q95_l = quantile(n_snps,probs = .025),
            n_snps_q95_u = quantile(n_snps,probs = .975))
