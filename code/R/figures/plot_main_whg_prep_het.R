library(tidyverse)
library(plyranges)
library(glue)
library(prismatic)
library(patchwork)
library(here)
library(data.table)
library(ggrastr)
source(here("code/R/project_defaults_shared.R"))

window_width <- 1e6
window_step <- 25e4

data_phenotype <- read_tsv(here("data/file_info.tsv")) |> 
  filter(!duplicated(sample_id)) |> 
  select(ind = sample_id, treatment) |> 
  mutate(treatment = replace_na(treatment, "mirleo"))

data_avg_ind <- read_tsv(here(glue("results/het/win_het_ind_all_w{window_width*1e-6}Mb_s{window_step*1e-3}kb.tsv.gz"))) |>
  filter(!is.na(avg_hom),
         start %% 1000000 == 0) |> 
  mutate(win_length = end - start + 1,
         het = 1 - avg_hom) |> 
  group_by(spec, ind) |> 
  summarise(avg_het_bp = sum(het * n_snps) / (sum(win_length)),
            avg_het_snp = sum(het * n_snps) / sum(n_snps)) |> 
  ungroup() |> 
  left_join(data_phenotype)

p1 <- data_avg_ind |> 
  ggplot(aes(x = spec, y = avg_het_bp)) +
  geom_jitter(height = 0, width = .33,
              shape = 21, color = clr_default[[1]],
              fill = clr_alpha(clr_default[[1]])) +
  labs(y = "Heterozygosity across samples (avg. over genome)",
       subtitle = "het by individuals (avg over genome)") +
  theme_ms() +
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.text.x = element_text(face = "italic"),
        plot.subtitle = element_text(color = "red"))

saveRDS(object = p1,
        here("results/img/R/p_het_ind_bp.Rds"))  

p1b <- data_avg_ind |> 
  ggplot(aes(x = spec_names[spec], y = avg_het_bp)) +
  geom_jitter(height = 0, width = .33,
              shape = 21, size = point_sz,
              aes( color = treatment,
              fill = after_scale(clr_alpha(color)))) +
  labs(y = "Heterozygosity") +
  scale_color_manual(values = c(clr_pheno, mirleo = clr_default[[2]]),
                     labels = c(lab_pheno, mirleo = "M. leonina"),
                     guide = "none") +
  coord_cartesian(ylim = c(0, .0018))+
  theme_ms() +
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.text.x = element_text(face = "italic"))

data_avg_ind |> 
  group_by(spec) |> 
  summarise(mean_het = mean(avg_het_bp),
            sd_het = sd(avg_het_bp)) |> 
  mutate(across(-spec, \(x){sprintf("%.6f", x)}))

saveRDS(object = p1b,
        here("results/img/R/p_het_ind_bp_pheno.Rds"))  

p2 <- data_avg_ind |> 
  ggplot(aes(x = spec, y = avg_het_snp)) +
  geom_jitter(height = 0, width = .33,
              size = point_sz,
              shape = 21, color = clr_default[[1]],
              fill = clr_alpha(clr_default[[1]])) +
  labs(y = "Heterozygosity across samples (avg. over SNPs)",
       subtitle = "het by individuals (avg over SNPs)") +
  theme_ms() +
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.text.x = element_text(face = "italic"),
        plot.subtitle = element_text(color = "red"))

saveRDS(object = p2,
        here("results/img/R/p_het_ind_snp.Rds"))  

data_avg_spec <- read_tsv(here(glue("results/het/win_het_by_spec_w{window_width*1e-6}Mb_s{window_step*1e-3}kb.tsv.gz")))

print_range <- \(x){
  rng <- range(x)
  str_c(c(sprintf("%.3f", rng[[1]])," – " ,sprintf("%.3f", rng[[2]])), collapse = "")
}

p3 <- data_avg_spec |> 
  ggplot(aes(x = spec_names[spec], y = win_het)) +
  geom_boxplot(color = clr_default[[1]],
               fill = clr_alpha(clr_default[[1]]),
               outlier.colour = "transparent") +
  geom_text(data = data_avg_spec |> 
              group_by(spec) |> 
              summarise(range = print_range(win_het)),
            aes(y = .003575, label = range),
            color = clr_default[[1]],
            family = fnt_sel,
            size = 4) +
  coord_cartesian(ylim = c(0, .0036)) +
  labs(y = "Heterozygosity across the genome (1Mb windows, avg. over samples)",
       subtitle = "het by bp, windows along genome (avg over samples)") +
  theme_ms()+
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.text.x = element_text(face = "italic"),
        plot.subtitle = element_text(color = "red"))

saveRDS(object = p3,
        here("results/img/R/p_het_whg_bp.Rds"))  

p4 <- data_avg_spec |> 
  ggplot(aes(x = spec_names[spec], y = avg_het)) +
  geom_boxplot(color = clr_default[[1]],
              fill = clr_alpha(clr_default[[1]]),
              outlier.colour = "transparent") +
  geom_text(data = data_avg_spec |> 
              group_by(spec) |> 
              summarise(range = print_range(avg_het)),
            aes(y = .75, label = range),
            color = clr_default[[1]],
            family = fnt_sel,
            size = 4) +
  coord_cartesian(ylim = c(0, .8))+
  labs(y = "Heterozygosity within SNPs (avg. over samples)",
       subtitle = "het by SNPs (avg over samples)") +
  theme_ms()+
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(family = fnt_sel),
        axis.text.x = element_text(face = "italic"),
        plot.subtitle = element_text(color = "red"))

saveRDS(object = p4,
        here("results/img/R/p_het_whg_snp.Rds"))  
