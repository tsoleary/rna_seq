# Script to load the data from the frontiers info and wrangle the data ---------

# load these required packages -------------------------------------------------
require(tidyverse)
source(here::here("functions.R"))


# load these files -------------------------------------------------------------

# differentially expressed genes results
setwd(here::here("cahan/results"))
res_cold <- read.csv(list.files(pattern = "cold_results"))
res_hot <- read.csv(list.files(pattern = "hot_results"))
deg <- full_join(res_cold, res_hot, by = "gene", suffix = c(".cold", ".hot"))

# phenotypic information and the gwas results
setwd(here::here("DRGP_GWAS"))
pheno <- read.csv("pheno.csv")
gwas_cold <- read.table("CTmin/gwas.top.4.annot", header = TRUE)
gwas_hot <- read.table("CTmax/gwas.top.4.annot", header = TRUE)
gwas_tb <- read.table("ThermalBreadth/gwas.top.4.annot", header = TRUE)

# the gwas with only 10^-5 to get effect sizes
gwas_cold_5 <- read.table("CTmin/gwas.top.annot", header = TRUE) %>%
  dplyr::filter(AvgMixedPval < 10^-5) %>%
  dplyr::select(ID, AvgEff) 
gwas_hot_5 <- read.table("CTmax/gwas.top.annot", header = TRUE) %>%
  dplyr::filter(AvgMixedPval < 10^-5) %>%
  dplyr::select(ID, AvgEff) 
gwas_tb_5 <- read.table("ThermalBreadth/gwas.top.annot", header = TRUE) %>%
  dplyr::filter(AvgMixedPval < 10^-5) %>%
  dplyr::select(ID, AvgEff) 

# transcription factor list 
dmel_tf <- read_delim(here::here("cahan/dmel_transcription_factors.txt"), 
                      delim = "\t")


# wrangle the data and annotation information ----------------------------------

# join the _5 to the _all to add the AvgEff in
gwas_cold <- gwas_cold %>%
  full_join(gwas_cold_5, by = "ID")
gwas_hot <- gwas_hot %>%
  full_join(gwas_hot_5, by = "ID")
gwas_tb <- gwas_tb %>%
  full_join(gwas_tb_5, by = "ID")

# annotated all the gene symbols and features from annot
gwas_cold <- gwas_dgrp_annot(gwas_cold)
gwas_hot <- gwas_dgrp_annot(gwas_hot)
gwas_tb <- gwas_dgrp_annot(gwas_tb)

# filter the gwas so there is only the pvalues that reach the cut off level
gwas_cold_5 <- gwas_cold %>%
  filter(AvgMixedPval < 10^-5)
gwas_hot_5 <- gwas_hot %>%
  filter(AvgMixedPval < 10^-5)
gwas_tb_5 <- gwas_tb %>%
  filter(AvgMixedPval < 10^-5)

# create tfs in hot and cold data
cold_tf <- res_cold %>%
  dplyr::filter(padj < 0.01) %>%
  dplyr::filter(gene %in% dmel_tf$SYMBOL) %>%
  dplyr::select(gene)

hot_tf <- res_hot %>%
  dplyr::filter(padj < 0.01) %>%
  dplyr::filter(gene %in% dmel_tf$SYMBOL) %>%
  dplyr::select(gene)


# integrate the two data sets --------------------------------------------------

# create a data.frame with info for the individual ma_plots 
deg_gwas_cold <- full_join(res_cold, gwas_cold, by = "gene")
deg_gwas_hot <- full_join(res_hot, gwas_hot, by = "gene")

# create a data.frame with info for the individual ma_plots 
deg_gwas_cold_5 <- full_join(res_cold, gwas_cold_5, by = "gene")
deg_gwas_hot_5 <- full_join(res_hot, gwas_hot_5, by = "gene")

# create a data.frame with all the info for the plotting later
total <- full_join(left_join(deg, 
                             gwas_cold, 
                             by = "gene", 
                             suffix = c("", ".cold")),
                   gwas_hot, 
                   by = "gene", 
                   suffix = c(".cold", ".hot"))

# group by the significance values
deg_grouped <- group_deg(total, pval_cut = 0.01, lfc_cut = 0)


# phenotype wrangling ----------------------------------------------------------
pheno <- pheno %>%
  pivot_longer(cols = starts_with("ct"), 
               names_to = "type",
               values_to = "temp") %>%
  separate(type, into = c("type", "sex") ,sep = "_")
