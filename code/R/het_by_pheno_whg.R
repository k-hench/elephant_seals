library(tidyverse)
library(glue)
library(prismatic)
library(patchwork)
library(here)
source(here("code/R/project_defaults.R"))

data_pheno <- read_tsv("results/pop/group_pheno_labeled.pop",
                       col_names = c("sample_id", "phenotype"),
                       col_types = "cc") |> 
  mutate(phenotype = factor(phenotype, levels = c("worms", "control", "mirleo")))

data_het <- map_dfr(c("mirang", "mirleo"), \(spec){read_tsv(here(glue("results/het/het_{spec}.tsv")), col_types = "cddid")}) |> 
  left_join(data_pheno, by = c(INDV = "sample_id"))

genome_size <- 2329128432

set.seed(42)
p1 <- data_het |> 
  ggplot(aes(x = phenotype, y = (N_SITES - `O(HOM)`)/N_SITES)) +
  labs(y = "Heterozygosity", subtitle = "Heterozygosity (fraction of SNPs)")

p2 <- data_het |> 
  ggplot(aes(x = phenotype, y = (N_SITES - `O(HOM)`)/genome_size)) +
  labs(y = "Heterozygosity", subtitle = "Heterozygosity (genome wide)") 

p3 <- data_het |> 
  group_by(phenotype) |>
  filter(row_number() == 1) |> 
  ungroup() |> 
  ggplot(aes(x = phenotype, y = N_SITES * 1e-6)) +
  labs(y = "n x 10^6", subtitle = "number of SNPs") 

p1 + p2 + p3 &
  theme_minimal(base_family = fnt_sel, base_size = 9) &
  geom_jitter(aes(color = phenotype, fill = after_scale(clr_alpha(color))),
              shape = 21, size = 2, height = 0, width = .2) &
  scale_color_manual(values = clr_pheno, guide = "none") &
  theme(axis.ticks = element_line(linewidth = .3),
        axis.line = element_line(linewidth = .3),
        axis.title.x = element_blank(),
        panel.grid = element_blank())

ggsave(filename = here(glue("results/img/het/het_genome_wide_by_phenotype.pdf")),
       # plot = p,
       width = 8, height = 3, device = cairo_pdf)

