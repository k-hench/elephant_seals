library(tidyverse)
library(prismatic)
library(here)
library(glue)
source(here("code/R/project_defaults_shared.R"))

read_pi <- \(part){read_csv(here(glue("results/pi/mirang_pi_dxy_{part}.tsv.gz")))}

data_pi <- str_pad(1:20, width = 2, pad = "0") |> 
  map_dfr(read_pi) |> 
  filter( (start - 1) %% 100000 == 0)

data_pi |> 
  summarise(mean_ang = mean(pi_mirang, na.rm = TRUE),
            sd_ang = sd(pi_mirang, na.rm = TRUE), 
            mean_leo = mean(pi_mirleo, na.rm = TRUE),
            sd_leo = sd(pi_mirleo, na.rm = TRUE)) |> 
  mutate(across(everything(), \(x){sprintf("%.4f", x)}))

avg_pi <- data_pi |> 
  select(sites, pi_mirang, pi_mirleo)|> 
  pivot_longer(cols = starts_with("pi"),
               names_to = "species",
               values_to = "pi",
               names_transform = \(str){str_remove(str,"pi_")}) |> 
  filter(complete.cases(pi)) |> 
  group_by(species) |> 
  summarise(genome_wide_avg_pi = sum(sites *pi) / sum(sites)) |> 
  ungroup()

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
  geom_point(data = avg_pi |> rename(pi = "genome_wide_avg_pi"),
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
  labs(y = "Nucleotide diversity (\U03C0)",
       x = "Species") +
  theme_ms() +
  theme(axis.text.x = element_text(face = "italic"))

saveRDS(object = p_pi,
        here("results/img/R/p_pi.Rds"))  

data_pi_plot |> 
  group_by(species) |> 
  summarise(min = min(pi, na.rm = TRUE),
            max = max(pi, na.rm = TRUE))
