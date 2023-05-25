# -----------------------------------------------------------
# script to sync {renv} with the conda environment `r_tidy`
# to create lock-file:
# conda activate r_tidy
# R
# renv::init()
# -----------------------------------------------------------
library("tidyverse")
library("here")
library("glue")
library("rlang")
library("readxl")
library("yaml")
library("lubridate")
library("patchwork")
library("cowplot")
library("prismatic")
library("ggnewscale")
library("ggforce")
library("geomtextpath")
library("ggstance")
library("ggridges")
library("sf")
library("ape")
library("tidytree")
library("renv")
library("BiocManager")
library("ggtree")


