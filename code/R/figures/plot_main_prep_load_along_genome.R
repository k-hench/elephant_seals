library(tidyverse)
library(patchwork)
library(glue)
library(prismatic)
library(here)
library(plyranges)

source(here("code/R/project_defaults_shared.R"))
specs <-  c("mirang", "mirleo")
spec_names <- c(mirang = "M. angustirostris", mirleo = "M. leonina")

import_genomes <- function(spec){
  read_tsv(here("data", "genomes", "filtered", glue("{spec}_filt.fa.gz.fai")),
           col_names = c("chr", "length", "offset", "linebases", "linewidth")) |> 
    mutate(spec = spec,
           idx = row_number(),
           end_pos = cumsum(length),
           start_pos = lag(end_pos, default = 0),
           eo = idx %% 2)
}

read_pop <- \(spec){read_tsv(here(glue("results/pop/inds_{spec}.pop")),
                             col_types = "c",
                             col_names = "sample") |>
    mutate(spec = spec)}
pops <- c("mirang", "mirleo") |> map_dfr(read_pop) 

genomes <- specs |> map_dfr(import_genomes)

g_starts <- dplyr::select(genomes, chr, ref = spec, start_pos) 

read_load <- \(sample, type){
  read_tsv(here(glue("results/mutation_load/snp_eff/by_ind/{type}/{sample}_{type}.bed.gz")),
           col_types = "ciic") |> 
    left_join(g_starts, by = c(`#CHROM` = "chr")) |> 
    mutate(gpos = FROM + start_pos) |> 
    select(chr = `#CHROM`, pos = FROM, gpos) |> 
    mutate(sample = sample, type = type)
}

samples <- dir("results/mutation_load/snp_eff/by_ind/expressed/") |> str_remove("_expressed.bed.gz")

data <- crossing(sample = samples,
                 type = c("expressed_anc", "masked", "fixed_anc")) |> 
  pmap_dfr(read_load) |> 
  left_join(pops)

load_labs <- c(masked = "inbreeding\nload",
               expressed = "segregating\nload", 
               fixed = "drift\nload")

data_counts <- data |> 
  group_by(chr, pos, gpos, type, spec) |> 
  count() |> 
  ungroup() |> 
  mutate(type = str_remove(type, "_anc"),
         load_label = factor(load_labs[type], levels = load_labs[c(2,1,3)]))

genome_mirang <- genomes |> 
  filter(spec == "mirang") |>
  mutate(mid_pos = (start_pos + end_pos) / 2)

clr_load_lab <- clr_load |> 
  set_names(nm = c(total = "total",
                   masked = load_labs[[1]],
                   fixed = load_labs[[3]],
                   expressed = load_labs[[2]])[names(clr_load)])

p <- ggplot() +
  geom_rect(data = genome_mirang , 
            aes(xmin = start_pos,
                xmax = end_pos, 
                ymin = -Inf,
                ymax = Inf, 
                fill = as.character(eo)),
            color = "transparent") +
  geom_point(data = data_counts,
             mapping =  aes(x = gpos,
                            y = n*2*(1.5- as.numeric(factor(spec))),
                            #color = load_label
                            color = spec
             ),
             size = 1.5, alpha = .5) +
  geom_hline(yintercept = 0, linewidth = .3, linetype = 3) +
  scale_x_continuous(name = "Genomic Position (Gb)",
                     labels = \(x){sprintf("%.1f", x * 1e-9)}#,
                     # sec.axis = sec_axis(name = "Scaffold Id",
                     #                     trans = identity,
                     #                     breaks = genome_mirang$mid_pos[1:17],
                     #                     label = genome_mirang$chr[1:17] |>
                     #                       str_remove("N[CW]_0723"))
                     )+
  scale_y_continuous("No. of Individuals with SNP",#"Sample count at SNP (n)",
                     labels = \(y){abs(y)}#,
                     # sec.axis = sec_axis(name = "Load Frequency",
                     #                     trans = identity,
                     #                     breaks = c(-20,-15,-10,-5,0,5,10,15,20),
                     #                     label = c(1,0.75, 0.5, 0.25, 0,
                     #                               0.25, 0.5, 0.75, 1))
  ) +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
  # scale_color_manual(values = clr_load_lab) +
  scale_color_manual("Species",
                     values = c("black", "gray50"),
                     labels = spec_names) +
  facet_grid(load_label ~ ., switch = "y") +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  coord_cartesian(ylim = c(-21, 21), expand = 0) +
  theme_ms() +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.text = element_text(face = "italic"))

ggsave(plot = p,
       here("results/img/load_along_genome_draft.pdf"),
       width = 12, height = 5, device = cairo_pdf)

saveRDS(object = p,
        here("results/img/R/p_load_along_genome.Rds"))  

data_genes <- plyranges::read_bed(here("results/genomes/mirang_genes.bed.gz"))

data_load_on_genes <- data_counts |>
  select(seqnames = chr, start = pos, end = pos, type, spec, n) |>
  group_by(spec) |> 
  nest() |> 
  mutate(granges = map(data, as_granges),
         overlaps = map(granges, \(gin){as_tibble(find_overlaps(gin, data_genes))})) |> 
  ungroup()

# data_load_on_genes |>
#   unnest(overlaps) |> 
#   group_by(type, spec, name) |> 
#   mutate(n_S = length(unique(end))) |> 
#   arrange(-n_S) |> View()

# data_genes_summary |>
#   select(type, spec, name, n_SNPs, n_SNPs_combined) |> 
#   pivot_wider(names_from = type, values_from = n_SNPs,values_fill = 0) |> View()
#   ggplot(aes(x = masked, y = expressed, color = fixed)) +
#   geom_jitter(alpha=.4, width = .3, height = .3) +
#   facet_wrap(spec ~ .) +
#   coord_equal()

data_genes_summary <- data_load_on_genes |>
  unnest(overlaps) |> 
  group_by(name, spec, seqnames) |> 
  mutate(n_alleles = if_else(type == "masked", n, 2*n),
         n_SNPs_combined = length(unique(end)),
         n_indv_combined = sum(n),
         n_alleles_combined = sum(n_alleles)) |> 
  ungroup() |> 
  group_by(type, spec, name) |> 
  summarise(n_SNPs_combined = n_SNPs_combined[[1]],
            n_indv_combined = n_indv_combined[[1]],
            n_alleles_combined = n_alleles_combined[[1]],
            n_SNPs = length(unique(end)),
            n_indv = sum(n),
            n_alleles = sum(n_alleles),
            spectrum = list(n)) |> 
  ungroup() |> 
  mutate(name = fct_reorder(name, -n_alleles_combined),
         load_label = factor(load_labs[type], levels = load_labs[c(2,1,3)]))

str_c("# ",
      c("Summary data for Elephant Seal mutation load distribution over genes within the genome annotation ('GCF_021288785.2_ASM2128878v3_genomic.gtf.gz')",
        "(header length: 12 lines)",
        "column details:",
        "species: Northern Elephant Seal (Mirounga angustirostris) / Southern Elephant Seal (Mirounga leonina)",
        "load_type: classification of load type followint the Vendrami et al. manuscript, depending on zygosity and fixation of the derived allele",
        "n_genes: number of genes for which at least one sample contains a SNP classified as the respective load_type",
        "n_SNPs_mean: average (mean) of the number of load SNPs/gene over all 20 samples",
        "n_SNPs_sd: standard deviation (SD) of the number of load SNPs/gene over all 20 samples",
        "n_indv_mean: average (mean) of the number of individual load scores SNPs/gene over all 20 samples (combination of n samples & n SNPs)",
        "n_indv_sd: standard deviation (SD) of the number of individual load scores SNPs/gene over all 20 samples (combination of n samples & n SNPs)",
        "n_alleles_mean: average (mean) of the number of alleles scored/gene over all 20 samples (this is equivalent to n_indv_* x 1 for inbreeding load, and n_indv_* x 2 for segregating and drift load)",
        "n_alleles_sd: standard deviation (SD) of the number of alleles scored/gene over all 20 samples (this is equivalent to n_indv_* x 1 for inbreeding load, and n_indv_* x 2 for segregating and drift load)")) |> 
  write_lines(here("results/tab/load_by_genes_summary.tsv"))


ggplot() +
  geom_histogram(data = data_genes_summary |> filter(!duplicated(name)),
                 aes(x = n_SNPs_combined),
                 binwidth = 1, color ="black", fill = clr_alpha("black"), linewidth = .3) +
  geom_histogram(data = data_genes_summary,
                 aes(x = n_SNPs_combined, y = after_stat(-count), color = type, fill = after_scale(clr_alpha(color))),
                 binwidth = 1, linewidth = .3) +
  facet_wrap(spec ~ .)

data_genes_summary |>
  group_by(spec) |> 
  filter(!duplicated(name)) |> 
  summarise(avg_n_snps_combined = mean(n_SNPs_combined),
            sd_n_snps_combined = sd(n_SNPs_combined))

data_genes_summary |>
  group_by(spec, type, load_label) |> 
  summarise(n_genes = length(n_SNPs),
            across(c(n_SNPs, n_indv, n_alleles),.fns=list(mean = mean, sd = sd),.names = "{.col}_{.fn}")) |> 
  mutate(load_label = str_replace(load_label, "\n", " "),
         species = spec_names[spec] |>  str_replace("M.", "Mirounga")) |> 
  ungroup() |> 
  select(species, load_type = load_label, n_genes:n_alleles_sd) |> 
  write_tsv(here("results/tab/load_by_genes_summary.tsv"),
            col_names = TRUE,
            append = TRUE)

str_c("# ",
      c("Mutation Load Tally by gene for Elephant Seals acording to the genome annotation ('GCF_021288785.2_ASM2128878v3_genomic.gtf.gz')",
        "(header length: 12 lines)",
        "column details:",
        "species: Northern Elephant Seal (Mirounga angustirostris) / Southern Elephant Seal (Mirounga leonina)",
        "load_type: classification of load type followint the Vendrami et al. manuscript, depending on zygosity and fixation of the derived allele",
        "gene_id: gene ID (field 'ID', extracted from 'GCF_021288785.2_ASM2128878v3_genomic.gtf.gz')",
        "n_SNPs_joint_load_types: number of SNPs on this gene when combining all three load types (which is different from the sum of the individual tallies, which would lead to double-counts of SNPs)",
        "n_indv_joint_load_types: number of individual load cases (combination of n samples & n SNPs) when combining all three load types (which is different from the sum of the individual tallies, which would lead to double-counts of SNPs)",
        "n_alleles_joint_load_types: number of load-alleles (n_indv_joint_load_types x 1 for inbreeding load, and n_indv_joint_load_types x 2 for segregating and drift load) when combining all three load types (which is different from the sum of the individual tallies, which would lead to double-counts of SNPs)",
        "n_SNPs: number of SNPs on with at least one samples scored for the respective load type at this gene",
        "n_indv: number of load cases scored for the respective load type at this gene (combination of n samples & n SNPs)",
        "n_alleles: number of alleles scored for the respective load type at this gene (n_indv x 1 for inbreeding load, and n_indv x 2 for segregating and drift load)")) |> 
  write_lines(here("results/tab/load_by_genes.tsv"))

data_genes_summary |> 
  mutate(load_label = str_replace(load_label, "\n", " "),
         species = spec_names[spec] |>  str_replace("M.", "Mirounga")) |> 
  select(-spectrum)|> 
  select(species,
         load_type = load_label,
         gene_id = name,
         n_SNPs_joint_load_types = n_SNPs_combined,
         n_indv_joint_load_types = n_indv_combined,
         n_alleles_joint_load_types = n_alleles_combined,
         n_SNPs,
         n_indv,
         n_alleles) |> 
  write_tsv(here("results/tab/load_by_genes.tsv"),
            col_names = TRUE,
            append = TRUE)
