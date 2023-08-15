library(tidyverse)
library(glue)
library(prismatic)
library(patchwork)

clr <- RColorBrewer::brewer.pal(3, "Set1") |> map(\(x){c(x, clr_desaturate(clr_lighten(x,.7),.1))}) |> unlist() |> color() |> 
  set_names(nm = str_c(rep(c("bot06", "bot10", "null"), each = 2), "-", rep(c("lgm", "nes"), 3)))

read_boostrap_sfs <- \(mod, bs_idx){
  bs_lab <- str_pad(bs_idx, width = 2, pad = 0)
  fl <- glue("results/demography/fastsimcoal/mirang_on_mirang/{mod}/bootstrap/bs_{bs_lab}/bestrun/mirang_on_mirang_{mod}_MAFpop0.txt")
  read_lines(fl, skip = 0, n_max = 2) |>
    str_replace_all("\t", " ") |>
    str_remove_all("d0_") |>
    str_split(" ") |>
    set_names(nm = c("allele_count", "snp_freq")) |>
    as_tibble() |>
    mutate(across(everything(), as.numeric)) |> 
    filter(!is.na(allele_count)) |> 
    pivot_wider(names_from = allele_count, values_from = "snp_freq") |> 
    mutate(mod = mod, bs_idx = bs_idx) |> 
    select(mod, bs_idx, everything())
}

data_exp_sfs <- cross_df(list(mod = c("bot06-lgm", "bot06-nes", "bot10-lgm", "bot10-nes", "null-lgm", "null-nes"),
              bs_idx = 1:50)) |> 
  arrange(mod, bs_idx) |> 
  pmap_dfr(read_boostrap_sfs)

p1 <- data_exp_sfs |> 
  pivot_longer(-c(mod,bs_idx), names_transform = as.numeric) |> 
  ggplot(aes(x = name, y = value, color = mod, group = str_c(mod, bs_idx))) +
  geom_step(linewidth = .2, alpha = .4) +
  facet_grid(mod ~ ., scales = "free") +
  labs(subtitle = "50 expected SFS (one per bootstrap run)") +
  scale_color_manual(values = clr) 

as.matrix(data_exp_sfs[,-c(1:2)])
# pca_data <- data_exp_sfs
pca_data <-  (data_exp_sfs |> filter(!grepl("null", mod)))
results <- prcomp(pca_data[,-c(1:3, 21:40)])

#reverse the signs
results$rotation <- -1*results$rotation

#display principal components

biplot(results, scale = 0)

loadings <- as_tibble(results$rotation) |> 
  mutate(allele_count = row_number())

p2 <- exp_var <- tibble(explained_var = c(results$sdev^2 / sum(results$sdev^2)),
       pc = seq_along(explained_var)) |> 
  ggplot(aes(x = pc, y = explained_var)) +
  geom_bar(stat = "identity", color = "gray30", fill = "gray90") +
  labs(subtitle = "PCA of the SFS")

scl <- .01
p3 <- as_tibble(cbind(pca_data[,c(1:2)], as_tibble(results$x))) |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = mod)) +
  geom_text(data = loadings, aes(x = PC1*scl, y = PC2*scl, label = allele_count)) +
  scale_color_manual(values = clr)

p1 + (p2 + p3 + plot_layout(ncol = 1, heights = c(.25, 1))) +
  plot_layout(guides = "collect") &
  theme_bw() &
  theme(legend.position = "bottom")

ggsave("~/Dropbox/David/elephant_seals/img/demography/SFS_PCA_no_null.pdf", width = 12, height = 8, device = cairo_pdf)
