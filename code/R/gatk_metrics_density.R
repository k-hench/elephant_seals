#!/usr/bin/env Rscript
# Rscript gatk_metrics_density.R <spec> <metrics.tsv>_
library(tidyverse)
library(glue)
library(EnvStats)
args <- commandArgs(trailingOnly = TRUE)
ref <- as.character(args[1])
tsv  <- as.character(args[2])

# reading in the metrics table
data <- read_tsv(tsv) |> mutate(log_FS = log10(FS))

# function to convert the raw observations
# to desnity table
density_n <- \(x, n = 301, adj = .6){
  rng <- range(x, na.rm = TRUE)
  if(any(is.infinite(range(rng) ))){rng <- c(-1,3)}  # catch log
  
  tibble(
      breaks = seq(rng[[1]], rng[[2]], length.out = n),
      density = demp(x = breaks,
                     obs = x,
                     density.arg.list = list(adjust = adj))
  )
}

# lookup the set of metric names
vars <- names(data)[-(1:2)]

# set the Name for the MultiQC report for each metric
nms <- c(MQ = "RMSMappingQuality (MQ)",
        QD = "QualByDepth (QD)",
        FS = "FisherStrand (FS)",
        log_FS = "FisherStrand (log10(FS))",
        SOR = " StrandOddsRatio (SOR)",
        MQRankSum = "MappingQualityRankSumTest (MQRankSum)",
        ReadPosRankSum = "ReadPosRankSumTest (ReadPosRankSum)")

# set the Description for the MultiQC report for each metric
desc <- c(MQ = "Root mean square mapping quality over all reads at site. Intended to include the standard deviation of the mapping qualities.",
        QD = "Variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.",
        FS = "Probability that there is strand bias at the site.",
        log_FS = "Probability that there is strand bias at the site. Here, the density of log10 of FS is displayed.",
        SOR = "Strand bias. SOR was created because FS tends to penalize variants that occur at the ends of exons which tend to only be covered by reads in one direction and FS gives those variants a bad score. SOR will take into account the ratios of reads that cover both alleles.",
        MQRankSum = "Compares the mapping qualities of the reads supporting the reference allele and the alternate allele (positive when mapping qualities alternate allele are higher than reference allele)",
        ReadPosRankSum = "Compares whether the positions of the reference and alternate alleles are different within the reads.")

# Export function which first parses the header for MultiQC
# and then adds the density table.
export_density <- \(var_in){
  var_in <- var_in
  proto_head <- c(as.character(glue("# parent_id: custom_section_{ref}")),
                  as.character(glue("# parent_name: 'raw SNP metrics ({ref})'")),
                  "# parent_description: 'Overview of the distribution of metrics for SNPs filtering'",
                  "# id: '{var_in}'",
                  "# section_name: '{nms[var_in]}'",
                  "# description: '{desc[var_in]}'",
                  '# pconfig: {{"height": 160}}')

  glue_env <- new.env()
  glue_env$proto_head <- proto_head
  glue_env$var_in <- var_in
  proto_head |>
    map_chr(glue, .envir = glue_env) |>
    write_lines(glue("../results/qc/snp_metrics/{ref}_{var_in}_dens_mqc.tsv"))

  data[[var_in]] |>
    density_n() |>
    set_names(nm = c(var_in, "density")) |>
    write_tsv(glue("../results/qc/snp_metrics/{ref}_{var_in}_dens_mqc.tsv"),
              append = TRUE,
              col_names = FALSE)
}

# apply the export function to all metrics
walk(vars, possibly(export_density, NULL))