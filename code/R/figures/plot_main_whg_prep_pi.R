library(tidyverse)
library(prismatic)
library(here)
library(glue)
source(here("code/R/project_defaults_shared.R"))

read_pi <- \(part){read_csv(here(glue("results/pi/mirang_pi_dxy_{part}.tsv.gz")))}

data_pi <- str_pad(1:20, width = 2, pad = "0") |> 
  map_dfr(read_pi)

avg_pi <- data_pi |> 
  summarise(mirang = sum(sites *pi_mirang) / sum(sites),
            mirleo = sum(sites *pi_mirleo) / sum(sites)) |> 
  pivot_longer(everything(),
               names_to = "species",
               values_to = "pi")

data_pi_plot <- data_pi |> 
  select(mirang = pi_mirang,
         mirleo = pi_mirleo) |> 
  pivot_longer(everything(),
               names_to = "species",
               values_to = "pi")

print_range <- \(x){
  rng <- range(x)
  str_c(c(rng[[1]]," – " ,sprintf("%.3f", rng[[2]])), collapse = "")
}

p_pi <- data_pi_plot |> 
  ggplot(aes(x = species, y = pi)) +
  geom_boxplot(outlier.color = NA,width = .4,
               color = clr_default[[2]],
               fill = clr_alpha(clr_default[[2]]),
               linewidth = .2) +
  geom_point(data = avg_pi,
             color = clr_default[[2]],
             fill = clr_alpha("white", .75),
             size = point_sz, shape  = 21) +
  # geom_text(data = data_pi_plot |> 
  #             group_by(species) |> 
  #             summarise(range = print_range(pi)),
  #           aes(y = .004, label = range),
  #           color = clr_default[[2]],
  #           family = fnt_sel,
  #           size = 3) +
  coord_cartesian(ylim = c(0,.004)) +
  scale_x_discrete(labels = \(x){spec_names[x]}) +
  labs(y = "Nucleotide diversity (\U03C0)") +
  theme_ms() +
  theme(axis.text.x = element_text(face = "italic"),
        axis.title.x = element_blank())

saveRDS(object = p_pi,
        here("results/img/R/p_pi.Rds"))  