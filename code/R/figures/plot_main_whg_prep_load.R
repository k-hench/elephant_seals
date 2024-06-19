library(tidyverse)
library(glue)
library(here)
library(prismatic)
library(patchwork)
source(here("code/R/project_defaults_shared.R"))

samples <- read_tsv(here("data/file_info.tsv")) |> 
  select(sample_id, spec, treatment) |> 
  filter(!duplicated(sample_id))

read_all_load <- \(load_type, sample_id){
  data <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}/{sample_id}_{load_type}.bed.gz")), col_types = "cddc")
  data_roh <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}_in_roh/{sample_id}_{load_type}_in_roh.bed.gz")),
                       col_names = c("#CHROM", "FROM", "TO"), col_types = "cdd")
  if(load_type == "masked"){
    data_anc <- data
    data_anc_roh <- data_roh
  } else {
    data_anc <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}_anc/{sample_id}_{load_type}_anc.bed.gz")), col_types = "cddc")
    data_anc_roh <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}_anc_in_roh/{sample_id}_{load_type}_anc_in_roh.bed.gz")),
                             col_names = c("#CHROM", "FROM", "TO"), col_types = "cdd")
  }
  
  tibble(sample_id = rep(sample_id, 4)) |> 
    left_join(samples) |> 
    mutate(load_type = rep(load_type, 4),
           snp_subset = c("all", "anc", "roh", "roh anc"),
           n_snps = c(nrow(data), nrow(data_anc), nrow(data_roh), nrow(data_anc_roh)),
           data = list(data, data_anc, data_roh, data_anc_roh))
}

load_types <- c("masked", "expressed", "fixed")
data_load <- expand_grid(sample_id = samples$sample_id,
                         load_type = load_types) |> 
  pmap_dfr(read_all_load) |> 
  mutate(load_type = factor(load_type, levels = load_types))

load_labs <- c("inbreeding\n ",
               "segregating", 
               "drift")

set.seed(42)
p1 <- data_load |> 
  filter(snp_subset == "anc") |> 
  mutate(load_label = c(masked = load_labs[[1]],
                        fixed = load_labs[[3]],
                        expressed = load_labs[[2]])[as.character(load_type)] |> 
           factor(levels = load_labs),
         load_nest = c(masked = "",
                       fixed = "realised load",
                       expressed = "realised load")[as.character(load_type)],
         treatment = factor(treatment, levels = names(clr_pheno)[c(1,2,4,3)])) |> 
  ggplot(aes(x = load_label, 
             y = n_snps,
             color = treatment)) +
  geom_jitter(shape = 21,
              height = 0, width = .25,
              size = point_sz,
              aes(fill = after_scale(clr_alpha(color))))+#, aes(color = sample_id == "160488")) +
  geomtextpath::geom_textsegment(inherit.aes = FALSE,
                                 data = tibble(y = -30, xmin = 1.7, xmax = 3.3),
                                 aes(y = y, yend = y, x = xmin, xend = xmax, label = "realised"),
                                 linewidth = .2, family = fnt_sel, size = 3) +
  facet_grid(. ~ spec_names[spec], scales = "free", switch = "y") +
  scale_color_manual("Phenotype",
                     values = clr_pheno,
                     labels = lab_pheno,
                     na.value = clr_default[2],
                     guide = guide_legend(title.position = "top")) +
  coord_cartesian(clip = FALSE,
                  ylim = c(0,235),
                  xlim = c(.4, 3.6),
                  expand = 0) +
  labs(y = "Load tally (no. SNPs)",
       subtitle = "Load Tally by Load Type",
       x = "Load Type") +
  theme_ms() +
  theme(#panel.background = element_rect(color = "gray80"),
    axis.line = element_line(),
    axis.ticks = element_line(color = "gray80"),
    # legend.margin = margin(5,0,5,0,"pt"),
    legend.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(face = "italic"),
    legend.position = "bottom",
    legend.title = element_text(hjust = .5))

saveRDS(object = p1,
        here("results/img/R/p_load_by_type_b.Rds"))

clr_lab <- c(clr_pheno, mirang = "black")
clr_load_lab <- clr_load |> 
  set_names(nm = c(total = "total",
                   masked = load_labs[[1]],
                   fixed = load_labs[[3]],
                   expressed = load_labs[[2]])[names(clr_load)])

p2 <- data_load |> 
  filter(snp_subset == "anc") |>
  mutate(sample_ord = factor(str_c(spec, replace_na(treatment, "mirang"), sample_id)),
         sample_lab = fct_reorder(glue("<span style='color:{clr_lab[replace_na(treatment, 'mirang')]}'>{sample_id}</span>"),
                                  as.numeric(sample_ord)),
         load_label = c(masked = load_labs[[1]],
                        fixed = load_labs[[3]],
                        expressed = load_labs[[2]])[as.character(load_type)] |> 
           factor(levels = load_labs)) |> 
  ggplot(aes(x = sample_lab, y = n_snps)) +
  geom_bar(stat = 'identity', aes(fill = load_label))+
  scale_fill_manual("Load type",
                    values = clr_load_lab,
                    labels = \(x){str_remove_all(str_remove(x, "load"),"\\n")}, 
                    guide = guide_legend(title.position = "top")) +
  coord_cartesian(#xlim = c(.4,30.6),
                  expand = 0) +
  labs(y = "Load tally (no. SNPs)",
       subtitle = "Individual Cummulative Load Tally",
       x = "Sample ID") +
  ggforce::facet_row(spec_names[spec]~., scales = "free_x",
                     space = "free") +
  theme_ms() +
  theme(axis.text.x = ggtext::element_markdown(angle = 90,
                                               hjust = 1,
                                               vjust = .5),
        legend.position = "bottom",
        legend.title = element_text(hjust = .5),
        strip.text = element_text(face = "italic"))

saveRDS(object = p2,
        here("results/img/R/p_load_by_ind_b.Rds"))  

data_load |> 
  filter(snp_subset == "anc") |>
  select(sample_id, spec, load_type, n_snps) |> 
  pivot_wider(names_from = load_type, values_from = n_snps) |> 
  mutate(total_load = masked + expressed + fixed,
         realized_load = expressed + fixed) |> 
  pivot_longer(-c(sample_id, spec), names_to = "load_type", values_to = "n_snps") |> 
  group_by(spec, load_type) |> 
  summarise(mean = mean(n_snps),
            sd = sd(n_snps)) |> 
  ungroup() |> 
  pivot_wider(names_from = spec, values_from = c(mean, sd)) |> 
  mutate(lab_ang = str_c(sprintf("%.1f", mean_mirang), " +/- ", sprintf("%.1f", sd_mirang)),
         lab_leo = str_c(sprintf("%.1f", mean_mirleo), " +/- ", sprintf("%.1f", sd_mirleo)))
