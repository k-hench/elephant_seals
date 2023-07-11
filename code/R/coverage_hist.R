library(tidyverse)
library(patchwork)
library(geomtextpath)
fnt_sel <- "Josefin sans"

read_cov_mask <- \(fl){
  read_tsv(fl,
           col_names = c("scaffold", "start", "end", "cov"),
           col_types = "ciii") |> 
    mutate(length = end - start) |> 
    group_by(cov) |> 
    summarise(basecount = sum(length)) |> 
    ungroup() |> 
    mutate(file = fl)
}

base_dir <- "results/qc/coverage/masks/"
fls <- c("105391_on_mirang_covmask.bed.gz",
         "507_on_mirang_covmask.bed.gz",
         "ES2829_on_mirang_covmask.bed.gz")

data <- str_c(base_dir, fls) |> 
  map_dfr(read_cov_mask)

file_info <- read_tsv("data/file_info.tsv") |> 
  select(sample_id, species, spec) |> 
  filter(!duplicated(sample_id))

data_c <- data |> 
  group_by( file ) |> 
  mutate(Gb_covered = cumsum(as.numeric(basecount)) * 1e-9) |> 
  ungroup() |> 
  mutate(file = str_remove(file, base_dir) |> str_remove(".bed.gz"),
         sample_id = str_remove(file, "_.*")) |> 
  left_join(file_info)

saveRDS(object = data_c, file = "~/Downloads/mirang_cov_hist.Rds")

p1 <- data_c |> 
  filter(cov < 75) |>
  # filter(cov < 2.5e3) |>
  ggplot(aes(x = cov, y = basecount * 1e-9)) +
  geom_bar(stat = "identity",
           aes(color = cov < 2, fill = after_scale(prismatic::clr_alpha(color))),
           linewidth = .15) +
  scale_color_manual(values = c(`TRUE` = "darkred", `FALSE` = "gray60"),
                     guide = "none") +
  facet_grid(spec + file ~ ., switch = "y") +
  labs(y = "covered refernce (Gb)", x = "ceverage (X)") +
  theme_minimal(base_family = fnt_sel) +
  theme(strip.placement = "outside")

mirang_total_length <- read_tsv("data/genomes/filtered/mirang_filt.fa.gz.fai",
         col_names = c("chr", "length", "offset", "linebases", "linewidth")) |> 
  pluck("length") |> 
  sum()

mirang_target = tibble(y =  rep(mirang_total_length * 1e-9, 2),
                       x = c(-Inf, Inf),
                       label = "total genome")

p2 <- data_c |> 
  filter(cov < 75) |>
  # filter(cov < 2.5e3) |>
  ggplot(aes(x = cov, y = Gb_covered)) +
  geom_textpath(data = mirang_target,
                aes(x = x, y= y, label = label),
                family = fnt_sel, color = "gray50",
                hjust = .13, linewidth = .2) +
  # geom_hline(yintercept =,
  #            linewidth = .3, color = "gray50") +
  geom_vline(xintercept = 1.5, linetype = 3, color = "darkred") +
  geom_area(color = "gray60", fill = prismatic::clr_alpha("gray60")) +
  facet_grid(spec + file ~ ., switch = "y") +
  labs(y = "cummulative genome cover (Gb)", x = "ceverage (X)") +
  theme_minimal(base_family = fnt_sel) +
  theme(strip.placement = "outside",
        strip.text = element_blank())

p1 + p2
ggsave("~/Downloads/geome_coverage_coverage_mask_2x.pdf", width = 12, height = 6.5, device = cairo_pdf)


data_c |> 
  group_by(file) |> 
  summarise(max_cov = max(cov))
