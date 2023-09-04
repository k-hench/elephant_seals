library(tidyverse)
library(glue)
library(here)
library(prismatic)
library(patchwork)
source(here("code/R/project_defaults_shared.R"))

samples <- read_tsv("data/file_info.tsv") |> 
  select(sample_id, spec, treatment) |> 
  filter(!duplicated(sample_id))

read_all_load <- \(load_type, sample_id){
  data <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}/{sample_id}_{load_type}.bed.gz")))
  data_roh <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}_in_roh/{sample_id}_{load_type}_in_roh.bed.gz")))
  if(load_type == "masked"){
    data_anc <- data
    data_anc_roh <- data_roh
  } else {
    data_anc <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}_anc/{sample_id}_{load_type}_anc.bed.gz")))
    data_anc_roh <- read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{load_type}_anc_in_roh/{sample_id}_{load_type}_anc_in_roh.bed.gz")))
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

# !! FIXED LOAD NEEDS TO BE CONSTANT WITHIN SPECIES 
# (the snps should be fixed after all)

data_load |> 
  ggplot(aes(x = load_type, y = n_snps, color = treatment)) +
  # ggbeeswarm::geom_beeswarm(shape = 21,
  #                           aes(fill = after_scale(clr_alpha(color))),
  #                           size = 2)+#, aes(color = sample_id == "160488")) +
  geom_jitter(shape = 21,
              height = 0, width = .25,
              aes(fill = after_scale(clr_alpha(color))),
              size = 2)+#, aes(color = sample_id == "160488")) +
  # geom_violin() +
  facet_grid(snp_subset ~ spec, scales = "free", switch = "y") +
  # ggrepel::geom_text_repel(mapping = aes(label = sample_id)) +
  scale_color_brewer(palette = "Set1", na.value = "black") +
  theme_minimal(base_family = fnt_sel) +
  theme(panel.background = element_rect(color = "gray80"),
        axis.line = element_line(),
        axis.ticks = element_line(color = "gray80"),
        strip.placement = "outside") +
  theme(legend.position = "bottom")

ggsave(here("results/img/load/load_snp_eff.pdf"),
       width = 8, height = 8, device = cairo_pdf)

set.seed(42)
p1 <- data_load |> 
  filter(snp_subset == "anc") |> 
  ggplot(aes(x = load_type, y = n_snps, color = treatment)) +
  geom_jitter(shape = 21,
              height = 0, width = .25,
              aes(fill = after_scale(clr_alpha(color))),
              size = 2)+#, aes(color = sample_id == "160488")) +
  facet_grid(. ~ spec, scales = "free", switch = "y") +
  scale_color_manual("Phenotype",
                     values = clr_pheno,
                     na.value = clr_default[2],
                     guide = guide_legend(title.position = "top")) +
  labs(y = "Load Tally (n SNPs)",
       subtitle = "Load Tally by Load Type") +
  theme_ms() +
  theme(#panel.background = element_rect(color = "gray80"),
        axis.line = element_line(),
        axis.ticks = element_line(color = "gray80"),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(hjust = .5))

clr_lab <- c(clr_pheno, mirang = clr_default[2])

p2 <- data_load |> 
  filter(snp_subset == "anc") |>
  mutate(sample_ord = factor(str_c(spec, replace_na(treatment, "mirang"), sample_id)),
         sample_lab = fct_reorder(glue("<span style='color:{clr_lab[replace_na(treatment, 'mirang')]}'>{sample_id}</span>"),
                                  as.numeric(sample_ord))) |> 
  ggplot(aes(x = sample_lab, y = n_snps)) +
  geom_bar(stat = 'identity', aes(fill = load_type))+
  scale_fill_manual("Load Type",
                    values = clr_load,
                     guide = guide_legend(title.position = "top")) +
  coord_cartesian(xlim = c(.4,30.6),
                  expand = 0) +
  labs(y = "Load Tally (n SNPs)",
       subtitle = "Individual Cummulative Load Tally") +
  theme_ms() +
  theme(axis.text.x = ggtext::element_markdown(angle = 90,
                                   hjust = 1,
                                   vjust = .5),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(hjust = .5))

p3 <- data_load |> 
  filter(snp_subset == "anc") |>
  mutate(sample_ord = factor(str_c(spec, replace_na(treatment, "mirang"), sample_id)),
         sample_lab = fct_reorder(glue("<span style='color:{clr_lab[replace_na(treatment, 'mirang')]}'>{sample_id}</span>"),
                                  as.numeric(sample_ord)),
         load_impact = if_else(load_type == "masked", n_snps, 2*n_snps)) |> 
  ggplot(aes(x = sample_lab, y = load_impact)) +
  geom_bar(stat = 'identity', aes(fill = load_type))+
  scale_fill_manual("Load Type",
                    values = clr_load,
                    guide = guide_legend(title.position = "top")) +
  coord_cartesian(xlim = c(.4,30.6),
                  expand = 0) +
  labs(y = "Load Score (n X scoring factor)",
       subtitle = "Load Score (scoring factor: Hom = 2, Het = 1)") +
  theme_ms() +
  theme(axis.text.x = ggtext::element_markdown(angle = 90,
                                   hjust = 1,
                                   vjust = .5),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(hjust = .5))

p_out <- p1 / (p2 + p3) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        plot.tag = element_text(family = fnt_sel))

ggsave(filename = here("results/img/load/load_snp_eff_anc_details.pdf"),
      plot = p_out,
      width = 9, height = 6, device = cairo_pdf)
