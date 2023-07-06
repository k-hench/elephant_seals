library(tidyverse)
library(here)
library(glue)

source(here("code/R/project_defaults.R"))

pop <- "mirang"
ref <- "mirang"

data_raw <- read_lines(here("results", "demography", "preview", glue("prev_{pop}_on_{ref}.txt")), skip = 11)

data <- tibble( pop = data_raw[1],
        values = str_split(data_raw[2], pattern = "\t")[[1]] |> str_remove_all("[\\(\\) ]")) |> 
  filter(values != "" ) |> 
  separate(values, into = c("n_haplotypes", "n_snps"), sep = ",", convert = TRUE)

p <- data |> 
  ggplot(aes(x = n_haplotypes, y = n_snps)) +
  geom_point() +
  labs(x = "n_haplotypes (2 x n)",
       title = glue("easySFS preview for {pop} on {ref}"),
       caption = glue("effect of down-sampling of samples on number of SNPs (to account for missing data)")) +
  theme_minimal(base_family = fnt_sel)

dir.create(here("results/img/demography"), showWarnings = FALSE)
ggsave(filename = here(glue("results/img/demography/preview_{pop}_on_{ref}.pdf")),
       plot = p,width = 5, height = 3.5, device = cairo_pdf)
