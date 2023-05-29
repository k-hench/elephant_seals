library(tidyverse)
library(patchwork)
library(here)

data_dir <- here("results", "genome_stats")
files <- str_c(data_dir, '/', dir(data_dir, pattern = ".tsv"))

data <- map_dfr(files,
                function(fl){
                  vroom::vroom(fl, col_types = "ccccdddD") |> 
                    mutate(spec = str_remove(fl, ".*/") |> str_remove(".tsv")) |> 
                    dplyr::select(spec, everything())}
                )

p0 <- data |> 
  ggplot(aes(y = spec, x = `Assembly Stats Scaffold N50`)) +
  geom_point() +
  scale_x_log10()

ggsave(filename = here("results", "img", "qc","genomes_n50.pdf"),
       plot = p0,
       height = 12,
       width = 6)

p1 <- data |> 
  mutate(spec = fct_reorder(spec, `Assembly Stats Total Sequence Length`)) |> 
  ggplot(aes(y = spec, x = `Assembly Stats Total Sequence Length`)) +
  geom_point()

ggsave(filename = here("results","img", "qc", "genomes_length.pdf"),
       plot = p1,
       height = 12,
       width = 6)

p2 <- data |> 
  mutate(`Assembly Stats Scaffold N50` = replace_na(`Assembly Stats Scaffold N50`, 0),
         spec = fct_reorder(spec,`Assembly Stats Scaffold N50`)) |> 
  ggplot(aes(y = spec, x = `Assembly Stats Scaffold N50`)) +
  geom_point() +
  scale_x_log10(limits = c(1, 1e9))

p3 <- data |>
  ggplot(aes(y = `Assembly Stats Total Sequence Length`,
             x = `Assembly Stats Scaffold N50`)) +
  geom_point() +
  scale_x_log10(limits = c(1, 1e9)) +
  scale_y_continuous("genome length", limits = c(2e9, 3.5e9))

p4 <- p2 / p3 + plot_layout(heights = c(1,.2))

ggsave(filename = here("results", "img", "qc", "genomes_sorted.pdf"),
       plot = p3,
       height = 12,
       width = 6)

write_tsv(data,
          file = here("results", "genome_stats.tsv"))