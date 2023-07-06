library(tidyverse)
library(here)
library(glue)
library(prismatic)

source(here("code/R/project_defaults.R"))

pop <- "mirang"
ref <- "mirang"

data_raw <- read_lines(here(glue("results/demography/sfs/{pop}_on_{ref}/fastsimcoal2/{pop}_MAFpop0.obs")))

data <- tibble(n_haplotypes = str_split(data_raw[2], "\t")[[1]] |> str_remove("d0_") |> as.numeric(),
               count = str_split(data_raw[3], " ")[[1]] |> as.numeric()) |> 
  mutate(maf = n_haplotypes / max(n_haplotypes))

p <- data |> 
  filter(maf <= .5) |> 
  ggplot(aes(x = maf, y = count)) +
  geom_bar(stat = "identity", color = "gray60", fill = clr_alpha("gray85"))+
  labs(title = glue(" folded SFS for {pop} on {ref}")) +
  theme_minimal(base_family = fnt_sel)

dir.create(here("results/img/demography"), showWarnings = FALSE)
ggsave(filename = here(glue("results/img/demography/sfs_{pop}_on_{ref}.pdf")),
       plot = p,width = 5, height = 3.5, device = cairo_pdf)
