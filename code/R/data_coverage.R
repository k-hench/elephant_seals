library(tidyverse)
library(here)
library(glue)

read_cov <- \(file){
  read_tsv(here(glue("results/qc/coverage/{file}"))) |> 
    rename( name  = "#rname") |> 
    mutate(w_depth = meandepth * endpos) |> 
    summarise(mean_depth = sum(w_depth) / sum(endpos),
              mean_cov = sum(covbases) / sum(endpos)) |> 
    mutate(sample = str_remove(file, "_on_mirang_coverage.tsv.gz"))
}


files <- dir("results/qc/coverage/", pattern = "gz")

data <- files |> 
  map_dfr(read_cov) |>
  mutate(spec = c("mirleo", "mirang")[grepl("ES", sample) + 1]) 


data |> 
  ggplot(aes(x = sample, y = mean_depth, fill = spec)) +
  geom_bar(stat= "identity")

data |> 
  group_by(spec) |> 
  summarise(across(starts_with("mean"),
                   .fns = list(mean = mean,
                               sd = sd,
                               median = median))) |> 
  write_tsv("~/Dropbox/Elephant seal paper/1 Data/wgs_coverage.tsv")
