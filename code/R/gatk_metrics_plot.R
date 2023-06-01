#!/usr/bin/env Rscript
# Rscript gatk_metrics_plot.R <spec>_
library(tidyverse)
library(patchwork)
args <- commandArgs(trailingOnly = TRUE)
ref <- as.character(args[1])

base_dir <- "../results/qc/snp_metrics/"
files <- dir(base_dir, pattern = ref)
files <- files[!grepl("_snp_metrics.tsv.gz", files)]

filter_values <- list(
  mirang = list(QD = 7.5, FS = 17.5, log_FS = log10(17.5), MQ = 55.0, SOR = 3.0, MQRankSum = c( -0.5, 0.5), ReadPosRankSum = c(-2.25, 2.25)),
  mirleo = list(QD = 4.0, FS = 10.0, log_FS = log10(10.0), MQ = 55.0, SOR = 3.0, MQRankSum = c(-0.5, 0.5), ReadPosRankSum = c(-2.5, 2.5)) )

filter_fun <- \(x, x_lab, ref){
  case_when(
    x_lab %in% c("FS", "log_FS", "SOR") ~ x < filter_values[[ref]][[x_lab]],
    x_lab %in% c("MQ", "QD") ~ x > filter_values[[ref]][[x_lab]],
    x_lab %in% c("MQRankSum", "ReadPosRankSum") ~ between(x, filter_values[[ref]][[x_lab]][1], filter_values[[ref]][[x_lab]][2])
  )
}

plt_dens <- \(idx){
  fl <- files[[idx]]
  x_lab <- str_remove(fl, "_dens_mqc.tsv") |> str_remove("^[a-z]{6}_" )
  ref <- str_sub(fl, 1, 6)
  dat <- read_tsv(str_c(base_dir, fl), skip = 7, col_names = c("x", "y"))
  dat |> 
    ggplot() +
    geom_area(aes(x = x,y = y),
              color = "gray70",
              fill = rgb(.7,.7,.7,.2)) +
    geom_area(data = dat |> 
                filter(filter_fun(x, x_lab, ref)),
              aes(x = x,y = y),
              color = "red",
              fill = rgb(.7,0,0,.2)) +
    geom_vline(xintercept = filter_values[[ref]][[x_lab]],
               color = "red", linetype = 3) +
    labs(x = x_lab)
}

p <- seq_along(files) |> 
  map(plt_dens) |> 
  wrap_plots(nrow = 2) +
  plot_annotation(title = str_c(str_sub(files[1],1,6), " SNP filter stats")) &
  theme_minimal() &
  theme(plot.title = element_text(hjust = .5),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 

ggsave(filename = str_c("../results/img/control/snp_metrics_", ref,".pdf"),
       plot = p,
       width = 8,
       height = 3,
       device = cairo_pdf)