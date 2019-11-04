# gwas dgrp --------------------------------------------------------------------

require(AnnotationDbi)
require(org.Dm.eg.db)
require(tidyverse)

# directories
gwas_directory <- "/Users/tsoleary/R/rna_seq/DRGP_GWAS/"
deg_directory <- "/Users/tsoleary/R/rna_seq/results/"

# import GWAS data
setwd(gwas_directory)
gwas_hot <- read.table("CTmax/gwas.top.annot", header = TRUE)
gwas_cold <- read.table("CTmin/gwas.top.annot", header = TRUE)

# import DEG data
setwd("/Users/tsoleary/R/rna_seq/results/")
deg_hot <- read.csv("Dm_DESeq2_hot_ctrl_results_with_normalized.csv")
deg_cold <- read.csv("Dm_DESeq2_cold_ctrl_results_with_normalized.csv")


# snps <- min_gwas_top$ID
# 
# # fly base requires the SNPs in a certain format
# chromosome <-  str_extract(snps, "[[:digit:]]?[[:alpha:]]")
# position <- str_replace(snps, "[[:digit:]]?[[:alpha:]]_", "") %>%
#   str_replace("_[[:alnum:]]+", "")
# 
# snps_flybase <- str_c(chromosome, 
#                       str_c(position, position, sep = "-"), sep = ":")

# get the FBgn# for all the snps
gwas_cold$FBgn <- str_extract(gwas_cold$GeneAnnotation, 
                                 "FBgn[[:digit:]]+")
gwas_hot$FBgn <- str_extract(gwas_hot$GeneAnnotation, 
                             "FBgn[[:digit:]]+")

# FBgn to gene_symbol
gwas_cold$gene <- as.character(mapIds(org.Dm.eg.db, 
                               keys = gwas_cold$FBgn, 
                               column = "SYMBOL", 
                               keytype = "FLYBASE",
                               multiVals = "first"))

gwas_hot$gene <- as.character(mapIds(org.Dm.eg.db, 
                                     keys = gwas_hot$FBgn, 
                                     column = "SYMBOL", 
                                     keytype = "FLYBASE",
                                     multiVals = "first"))

comb_cold <- full_join(deg_cold, gwas_cold, by = "gene")

comb_hot <- full_join(deg_hot, gwas_hot, by = "gene")

# get the median or average for each treatment then do the 
cold <- comb_cold %>% 
  dplyr::select(contains("_"), gene, ID) %>%
  pivot_longer(contains("_"), names_to = "group", values_to = "expression") %>%
  mutate(group = str_replace(group, "_[[:digit:]]*$", "")) %>%
  group_by(gene, group, ID) %>%
  summarize(expression = median(expression, na.rm = TRUE)) %>%
  filter(expression > 0) %>%
  pivot_wider(names_from = group, values_from = expression) %>%

ggplot(cold) + 
  geom_point(aes(x = con, y = cold)) + 
  scale_x_log10() +
  scale_y_log10()
  

         