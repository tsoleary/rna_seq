# ------------------------------------------------------------------------------
# Dan - Huang et al 2015 25°C exp and CTmin CTmax correlations
# March 16, 2021
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)
require(broom)
require(here)

# Source functions
source(here::here("dan/functions.R"))

# Load data --------------------------------------------------------------------

# Load frontiers phenotype data
pheno <- read_csv(here::here("dan/pheno.csv"))

# Load Huang et al 2015 normalized expression data
exp_f <- read_delim(here::here("dan/dgrp.array.exp.female.txt"), 
                    delim = " ") 
colnames(exp_f)[2:length(colnames(exp_f))] <- 
  paste0(colnames(exp_f)[2:length(colnames(exp_f))], ":female")

exp_m <- read_delim(here::here("dan/dgrp.array.exp.male.txt"), 
                    delim = " ")
colnames(exp_m)[2:length(colnames(exp_m))] <- 
  paste0(colnames(exp_m)[2:length(colnames(exp_m))], ":male")

# Join the male and female expression data together
exp <- full_join(exp_f, exp_m, by = "gene") 

# Wrangle data -----------------------------------------------------------------

# Keep only expression data that has phenotype data
exp <- exp %>% 
  select(gene, contains(as.character(pheno$line_DGRP)))

# DGRP lines that have expression data
exp_lines <- str_remove_all(str_remove_all(colnames(exp), "line_"), 
                            ":[0-9]:[a-z]+")

# Keep only phenotype data that has expression data
metadata <- pheno %>%
  filter(line_DGRP %in% exp_lines)

# Play with data ---------------------------------------------------------------


# Get exp data for a specific gene with CTmin/max values for males & females
df <- exp_data_tidy("FBgn0000014", pheno)


# Run a linear model on that gene for CTmin/max values for males & females
df %>%
  group_by(sex_exp, limit) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(norm_exp ~ temp, data = .x)),
    tidied = map(fit, tidy),
  ) %>%
  unnest(tidied)


# Plot the normalized expression against CTmin and CTmax for a gene
plot_min_max_gene("FBgn0000014")


# Plot a bunch of the data -----------------------------------------------------

plot_list <- vector(mode = "list", length = nrow(exp))

for (i in 1:nrow(exp)) {
  plot_list[[i]] <- plot_min_max_gene(exp$gene[i])
}

pdf("all_plots.pdf", width = 10, height = 6)

for(i in 1:nrow(exp)){
  print(plot_list[[i]])
}

dev.off() # pdf file should appear in working directory



