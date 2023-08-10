library(tidyverse)
library(plyranges)
library(glue)
library(prismatic)
library(patchwork)
library(here)
library(data.table)
library(ggrastr)
source(here("code/R/project_defaults.R"))

import_genomes <- function(spec){
  read_tsv(here("data", "genomes", "filtered", glue("{spec}_filt.fa.gz.fai")),
           col_names = c("chr", "length", "offset", "linebases", "linewidth")) |> 
    mutate(spec = spec,
           idx = row_number(),
           end_pos = cumsum(length),
           start_pos = lag(end_pos, default = 0),
           eo = idx %% 2)
}

genome_tib <- import_genomes("mirang")
genome <- genome_tib  |> 
  mutate(start = 0) |> 
  dplyr::select(seqnames = chr, start, end = length) |> 
  as_granges()

window_width <- 1e6
window_step <- 25e4

# window_width <- 2.5e5
# window_step <- 5e4
window_mirang <- slide_ranges(genome, width = window_width, step = window_step)

convert_012_to_hom <- \(mat){abs(1 - mat)}
convert_012_to_het <- \(mat){1 - abs(1 - mat)}

avg_het_by_snp <- \(spec#, n_total = 20
                    ){
  # set.seed(42)
  # sampl_sel <- sample(1:n_total, size = 10, replace = FALSE)
  geno_mat <- (fread(glue("results/genotyping/012/mirang_filtered_{spec}_012.tsv.gz"), sep = "\t"))[,#sampl_sel,
                                                                                                    -1] |> as.matrix()
  geno_mat[geno_mat == -1] <- NA
  cum_het <- colSums(convert_012_to_het(geno_mat), na.rm = TRUE)
  n_typed <- colSums(!is.na(geno_mat), na.rm = TRUE)
  avg_het <- cum_het / n_typed
  avg_het
}

avg_het <- list(
  mirang = avg_het_by_snp("mirang"),
  mirleo = avg_het_by_snp("mirleo")
)

win_avg_het <- \(spec){
  snp_ids <- read_tsv(glue("results/genotyping/012/mirang_filtered_{spec}_012.012.pos"),
                      col_names = c("chr", "pos")) |>
    mutate(snp = str_c(chr,"_",pos))
  
  granges_012 <- snp_ids |> 
    mutate(mean_het = avg_het[[spec]],
           seqnames = str_remove(snp, "_[0-9]*$"),
           start = str_remove_all(snp, ".*_") |> as.integer()) |> 
    mutate(end = start) |> 
    select(-c(chr, snp)) |> 
    as_granges()
  
  join_overlap_inner(window_mirang, granges_012) |>
    as_tibble() |>
    group_by(seqnames, start, end, strand, partition) |>
    summarise(avg_het = mean(mean_het, na.rm = TRUE),
              win_het = sum(mean_het, na.rm = TRUE) / (end[1] - start[1])) |>
    mutate(mid = (start + end + 2) /2) |>
    left_join(genome_tib |>
                dplyr::select(seqnames = chr, gstart = start_pos, eo)) |>
    mutate(gpos = mid + gstart,
           spec = spec) |>
    ungroup()
}

windowed_avg_het <- c("mirang", "mirleo") |> 
  map_dfr(win_avg_het)

# windowed_avg_het |>
#   ggplot() +
#   geom_rect(data = genome_tib |> select(-spec),
#             aes(xmin = start_pos,
#                 xmax = end_pos,
#                 ymin = -Inf,
#                 ymax = Inf,
#                 fill = as.character(eo)),
#             color = "transparent") +
#   geom_line(aes(x = gpos, y = avg_het, group = str_c(seqnames, spec), color = spec), linewidth = .25, alpha = .2) +
#   # geom_point(aes(x = gpos, y = avg_het, group = seqnames, color = spec), size = .15, shape = 19, alpha = .2) +
#   scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
#   scale_x_continuous(expand = c(0, 0)) +
#   # facet_grid(spec ~ .) +
#   # coord_cartesian(xlim = c(0, max(genome_tib$end_pos)),
#   #                 ylim = c(-.1 , 1.1), 
#   #                 expand = 0) +
#   scale_color_manual(values = clrs, guide = "none") +
#   labs(subtitle = glue("heterozygosity averaged by species (overlay)")) +
#   theme_minimal(base_family = fnt_sel, base_size = 9) +
#   theme(axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         # axis.text.y = element_blank(),
#         panel.grid = element_blank())
# 
# ggsave(filename = here(glue("results/img/het/win_avg_het_overlay_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.png")),
#        # plot = p,
#        width = 10, height = 3, bg = "white")

p1 <- windowed_avg_het |> 
  filter(seqnames %in% genome_tib$chr[2:17]) |> 
  ggplot(aes(x = spec, y = avg_het, color = spec)) +
  geom_boxplot(aes(fill = after_scale(clr_alpha(color)))) +
  facet_wrap(seqnames ~ ., nrow = 1) +
  scale_color_manual(values = clrs, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(subtitle = "by autosome") +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank())

p2 <- windowed_avg_het |> 
  filter(seqnames %in% genome_tib$chr[2:17]) |> 
  ggplot(aes(x = spec, y = avg_het, color = spec)) +
  geom_boxplot(aes(fill = after_scale(clr_alpha(color)))) +
  scale_color_manual(values = clrs, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "heterozygosity", subtitle = glue("combined autosomes")) +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank())

p2 + p1 +
  plot_layout(widths = c(.1, 1)) +
  plot_annotation(subtitle = glue("heterozygosity within SNPs only (w{sprintf('%.0fMb', window_width /1e6)}, s{sprintf('%.0fkb', window_step /1e3)})"),
                  theme = theme(title = element_text(family = fnt_sel)))

ggsave(filename = here(glue("results/img/het/win_avg_het_autosomes_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}_10inds.png")),
       # plot = p,
       width = 15, height = 4, bg = "white")

p3 <- windowed_avg_het |> 
  filter(seqnames %in% genome_tib$chr[2:17]) |> 
  ggplot(aes(x = spec, y = win_het, color = spec)) +
  geom_boxplot(aes(fill = after_scale(clr_alpha(color)))) +
  facet_wrap(seqnames ~ ., nrow = 1) +
  scale_color_manual(values = clrs, guide = "none") +
  labs(subtitle = "by autosome") +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank())

p4 <- windowed_avg_het |> 
  filter(seqnames %in% genome_tib$chr[2:17]) |> 
  ggplot(aes(x = spec, y = win_het, color = spec)) +
  geom_boxplot(aes(fill = after_scale(clr_alpha(color)))) +
  scale_color_manual(values = clrs, guide = "none") +
  labs(y = "heterozygosity", subtitle = glue("combined autosomes")) +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank())

p4 + p3 +
  plot_layout(widths = c(.1, 1)) +
  plot_annotation(subtitle = glue("heterozygosity averaged over genome size (w{sprintf('%.0fMb', window_width /1e6)}, s{sprintf('%.0fkb', window_step /1e3)})"),
                  theme = theme(title = element_text(family = fnt_sel))) &
  coord_cartesian(ylim = c(0, .0062)) 

ggsave(filename = here(glue("results/img/het/win_avg_het_by_bp_autosomes_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}_10inds.png")),
       # plot = p,
       width = 15, height = 4, bg = "white")


# run_window_smlh <- \(spec = "mirang", smpl_idx){
#   snp_ids <- read_tsv(glue("results/genotyping/012/mirang_filtered_{spec}_012.012.pos"),
#                       col_names = c("chr", "pos")) |>
#     mutate(snp = str_c(chr,"_",pos))
#   inds <- read_lines(glue("results/genotyping/012/mirang_filtered_{spec}_012.012.indv"))
#   
#   geno_mat <- (fread(glue("results/genotyping/012/mirang_filtered_{spec}_012.tsv.gz"), sep = "\t")[smpl_idx,]) |> as.matrix()
#   geno_mat[geno_mat == -1] <- NA
#   
#   granges_012 <- snp_ids |> 
#     mutate(hom = convert_012_to_hom(geno_mat)[1,-1],
#            mean_het = avg_het[[spec]],
#            seqnames = str_remove(snp, "_[0-9]*$"),
#            start = str_remove_all(snp, ".*_") |> as.integer()) |> 
#     mutate(end = start) |> 
#     select(-c(chr, snp)) |> 
#     as_granges()
# 
#   join_overlap_inner(window_mirang, granges_012) |>
#     as_tibble() |>
#     group_by(seqnames, start, end, strand, partition) |>
#     summarise(avg_hom = mean(hom, na.rm = TRUE),
#               n_snps = n(),
#               sMLH = (1 - mean(hom, na.rm = TRUE)) / mean(mean_het, na.rm = TRUE)) |>
#     mutate(mid = (start + end + 2) /2) |>
#     left_join(genome_tib |>
#                 dplyr::select(seqnames = chr, gstart = start_pos, eo)) |>
#     mutate(gpos = mid + gstart,
#            ind = inds[smpl_idx],
#            spec = spec) |>
#     ungroup()
#   # avg_hom
# }
# 
# data_smlh <- run_window_smlh("mirang", 1)
# 
# data_smlh |>
#     ggplot() +
#     geom_rect(data = genome_tib,
#               aes(xmin = start_pos,
#                   xmax = end_pos,
#                   ymin = -Inf,
#                   ymax = Inf,
#                   fill = as.character(eo)),
#               color = "transparent") +
#     geom_line(aes(x = gpos, y = sMLH, group = seqnames, color = spec), linewidth = .25) +
#     geom_point(aes(x = gpos, y = sMLH, group = seqnames, color = spec), size = .15, shape = 19, alpha = .2) +
#     scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
#     scale_x_continuous(expand = c(0, 0)) +
#     # coord_cartesian(xlim = c(0, max(genome_tib$end_pos)),
#     #                 ylim = c(-.1 , 1.1), 
#     #                 expand = 0) +
#     scale_color_manual(values = clrs, guide = "none") +
#     # labs(subtitle = glue("{inds[smpl_idx]} ({spec})")) +
#     theme_minimal(base_family = fnt_sel, base_size = 9) +
#     theme(axis.title.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           # axis.text.y = element_blank(),
#           panel.grid = element_blank())
# 
# all_smlh_windows <- tibble(spec = rep(c("mirang", "mirleo"), c(20, 10)),
#                       smpl_idx = c(1:20, 1:10)) |> 
#   pmap_dfr(run_window_smlh)
# 
# 
# all_smlh_windows |> 
#   ggplot() +
#   geom_rect(data = genome_tib, 
#             aes(xmin = start_pos,
#                 xmax = end_pos, 
#                 ymin = -Inf,
#                 ymax = Inf, 
#                 fill = as.character(eo)),
#             color = "transparent") +
#   geom_hline(yintercept = 0:1, color = rgb(0,0,0,.1), linewidth = .2) +
#   geom_line(aes(x = gpos, y = sMLH, group = seqnames, color = spec), linewidth = .25)  +
#   geom_point(aes(x = gpos, y = sMLH, group = seqnames, color = spec), size = .15, shape = 19, alpha = .2) +
#   # scale_y_continuous(breaks = 0:1) +
#   scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') + 
#   coord_cartesian(xlim = c(0, max(genome_tib$end_pos)),
#                   ylim = c(-.1 , 5.1),
#                   expand = 0) +
#   scale_color_manual(values = clrs, guide = "none") + 
#   # labs(subtitle = glue("{inds[smpl_idx]} ({spec})")) +
#   facet_grid(ind ~ ., switch = "y") +
#   theme_minimal(base_family = fnt_sel, base_size = 9) +
#   theme(#axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         # axis.text.y = element_blank(),
#         panel.grid = element_blank(),
#         strip.placement = "outside",
#         strip.text.y.left = element_text(angle = 0))
# 
# ggsave(filename = here(glue("results/img/het/win_smlh_ind_all_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.png")),
#        # plot = p,
#        width = 10, height = 15, bg = "white")
# 
# write_tsv(all_smlh_windows,
#           here(glue("results/het/win_smlh_ind_all_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.tsv.gz")))


# consider phenotype within mirang ========================================================================
clr_pheno <-  RColorBrewer::brewer.pal(3, "Set1") |> 
  set_names(nm = c("worms", "control", "mirleo"))
data_pheno <- read_tsv("results/pop/group_pheno_labeled.pop",
                       col_names = c("sample_id", "phenotype"),
                       col_types = "cc") |> 
  mutate(phenotype = factor(phenotype, levels = c("worms", "control", "mirleo")))


avg_het_by_pheno <- \(spec, pheno){
  ind_idx <- tibble(sample_id = read_lines(glue("results/genotyping/012/mirang_filtered_{spec}_012.012.indv")),
                 idx = seq_along(sample_id)) |> 
    left_join(data_pheno) |> 
    filter(phenotype == pheno) |> 
    pluck("idx")         
  
  # set.seed(42)
  # sampl_sel <- sample(1:n_total, size = 10, replace = FALSE)
  geno_mat <- (fread(glue("results/genotyping/012/mirang_filtered_{spec}_012.tsv.gz"), sep = "\t"))[ind_idx, -1] |> as.matrix()
  geno_mat[geno_mat == -1] <- NA
  cum_het <- colSums(convert_012_to_het(geno_mat), na.rm = TRUE)
  n_typed <- colSums(!is.na(geno_mat), na.rm = TRUE)
  avg_het <- cum_het / n_typed
  avg_het
}

# read_het <- \(spec){
#   read_tsv(here(glue("results/het/het_{spec}.tsv")), col_types = "cddid") |> 
#     left_join(data_pheno, by = c(INDV = "sample_id"))
#   }
# 
# data_het_global <- c("mirang", "mirleo") |> 
#   map_dfr(read_het)
# 
# data_het_global |> 
#   ggplot(aes(x = phenotype, y = `O(HOM)`/N_SITES, color = phenotype)) +
#   geom_boxplot(aes(fill = after_scale(clr_alpha(color)))) +
#   scale_color_manual(values = clr_pheno, guide = "none") +
#   theme_minimal(base_family = fnt_sel, base_size = 9) +
#   theme(axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         # axis.text.y = element_blank(),
#         panel.grid = element_blank())

avg_het_pheno <- list(
  control = avg_het_by_pheno("mirang", pheno = "control"),
  worms = avg_het_by_pheno("mirang", pheno = "worms"),
  mirleo = avg_het_by_pheno("mirleo", pheno = "mirleo")
)

win_avg_het_pheno <- \(spec, pheno){
  snp_ids <- read_tsv(glue("results/genotyping/012/mirang_filtered_{spec}_012.012.pos"),
                      col_names = c("chr", "pos")) |>
    mutate(snp = str_c(chr,"_",pos))
  
  granges_012 <- snp_ids |> 
    mutate(mean_het = avg_het_pheno[[pheno]],
           seqnames = str_remove(snp, "_[0-9]*$"),
           start = str_remove_all(snp, ".*_") |> as.integer()) |> 
    mutate(end = start) |> 
    select(-c(chr, snp)) |> 
    as_granges()
  
  join_overlap_inner(window_mirang, granges_012) |>
    as_tibble() |>
    group_by(seqnames, start, end, strand, partition) |>
    summarise(avg_het = mean(mean_het, na.rm = TRUE),
              win_het = sum(mean_het, na.rm = TRUE) / (end[1] - start[1])) |>
    mutate(mid = (start + end + 2) /2) |>
    left_join(genome_tib |>
                dplyr::select(seqnames = chr, gstart = start_pos, eo)) |>
    mutate(gpos = mid + gstart,
           spec = spec,
           pheno = pheno) |>
    ungroup()
}

windowed_avg_het_pheno <- map2_dfr(c("mirang", "mirang", "mirleo"),
                                   c("control", "worms", "mirleo"),
                                   win_avg_het_pheno) |> 
  mutate(pheno = factor(pheno, levels = c("worms", "control", "mirleo")))


p5 <- windowed_avg_het_pheno |> 
  filter(seqnames %in% genome_tib$chr[2:17]) |> 
  ggplot(aes(x = pheno, y = avg_het, color = pheno)) +
  geom_boxplot(aes(fill = after_scale(clr_alpha(color)))) +
  facet_wrap(seqnames ~ ., nrow = 1) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1"), guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(subtitle = glue("heterozygosity in autosomes (w{sprintf('%.0fMb', window_width /1e6)}, s{sprintf('%.0fkb', window_step /1e3)})")) +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank())

p6 <- windowed_avg_het_pheno |> 
  filter(seqnames %in% genome_tib$chr[2:17]) |> 
  ggplot(aes(x = pheno, y = avg_het, color = pheno)) +
  geom_boxplot(aes(fill = after_scale(clr_alpha(color)))) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1"), guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(subtitle = glue("ALL autosomes")) +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank())

p6 + p5 + plot_layout(widths = c(.1, 1))

ggsave(filename = here(glue("results/img/het/win_avg_het_autosomes_pheno_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}_10inds.png")),
       # plot = p,
       width = 14, height = 4, bg = "white")

windowed_avg_het_pheno |>
  ggplot() +
  geom_rect(data = genome_tib |> select(-spec),
            aes(xmin = start_pos,
                xmax = end_pos,
                ymin = -Inf,
                ymax = Inf,
                fill = as.character(eo)),
            color = "transparent") +
  geom_line(aes(x = gpos, y = avg_het, group = str_c(seqnames, pheno), color = pheno), linewidth = .25, alpha = .2) +
  geom_point(aes(x = gpos, y = avg_het, group = seqnames, color = pheno), size = .15, shape = 19, alpha = .2) +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(pheno ~ ., switch = "y") +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1"), guide = "none") +
  labs(subtitle = glue("heterozygosity averaged by phenotype (based on SNPs only)"),
       y = "heterozygosity") +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.placement = "outside",
        # axis.text.y = element_blank(),
        panel.grid = element_blank())

ggsave(filename = here(glue("results/img/het/win_avg_het_pheno_snps_only_w{sprintf('%.0fMb', window_width /1e6)}_s{sprintf('%.0fkb', window_step /1e3)}.png")),
       # plot = p,
       width = 10, height = 8, bg = "white")

data_het_ind <- read_tsv(here(glue("results/het/win_het_ind_all_w1Mb_s250kb.tsv.gz"))) |> 
  left_join(data_pheno, by = c(ind = "sample_id"))

data_het_ind |> 
  rename(phenotype = "pheno") |> 
  filter(seqnames == genome_tib$chr[7]) |> 
  ggplot() +
  geom_rect(data = tibble(xmin = 8.125e7, xmax = 8.35e7),
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax =Inf),
            fill = "gray90") + 
  # geom_rect(data = genome_tib |> select(-spec),
  #           aes(xmin = start_pos,
  #               xmax = end_pos,
  #               ymin = -Inf,
  #               ymax = Inf,
  #               fill = as.character(eo)),
  #           color = "transparent") +
  geom_line(aes(x = gpos - gstart, y = 1 - avg_hom, group = str_c(seqnames, pheno, ind), color = pheno), linewidth = .25, alpha = .2) +
  geom_point(aes(x = gpos - gstart, y = 1 - avg_hom, group = seqnames, color = pheno), size = .15, shape = 19, alpha = .2) +
  scale_fill_manual(values = c(`0` = rgb(0,0,0,.1), `1` = rgb(0,0,0,0)), guide = 'none') +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(pheno ~ .) +
  # coord_cartesian(xlim = c(0, max(genome_tib$end_pos)),
  #                 ylim = c(-.1 , 1.1), 
  #                 expand = 0) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1"), guide = "none") +
  labs(subtitle = glue("heterozygosity averaged by species (overlay)")) +
  theme_minimal(base_family = fnt_sel, base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        panel.grid = element_blank())