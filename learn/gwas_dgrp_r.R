# gwas dgrp --------------------------------------------------------------------

require(AnnotationDbi)
require(org.Dm.eg.db)
require(tidyverse)

# directories
gwas_directory <- "/Users/tsoleary/R/rna_seq/DRGP_GWAS/"
deg_directory <- "/Users/tsoleary/R/rna_seq/results/"

# import GWAS data
setwd(gwas_directory)
gwas<- read.table("CTmin/gwas.top.annot", header = TRUE)

# import DEG data
setwd("/Users/tsoleary/R/rna_seq/results/")
deg <- read.csv("Dm_DESeq2_cold_ctrl_results_with_normalized.csv")

# get the FBgn# for all the snps
gwas$FBgn <- str_extract(gwas$GeneAnnotation, 
                         "FBgn[[:digit:]]+") 

# FBgn to gene_symbol
gwas$gene <- as.character(mapIds(org.Dm.eg.db, 
                          keys = gwas$FBgn, 
                          column = "SYMBOL", 
                          keytype = "FLYBASE",
                          multiVals = "first"))

comb <- full_join(deg, gwas, by = "gene")

# get the median or average for each treatment then do the 
comb_avg <- comb %>% 
  dplyr::select(contains("_"), gene, ID, padj) %>%
  pivot_longer(contains("_"), names_to = "group", values_to = "expression") %>%
  mutate(group = str_replace(group, "_[[:digit:]]*$", "")) %>%
  group_by(gene, group, ID, padj) %>%
  summarize(expression = median(expression, na.rm = TRUE)) %>%
  filter(expression > 0) %>%
  pivot_wider(names_from = group, values_from = expression)

comb_sort <- comb_avg %>%
  mutate(g = case_when(is.na(ID) & padj < 0.05 ~ "DEG",
                       is.na(ID) & padj >= 0.05 ~ "all",
                       is.na(ID) & is.na(padj) ~ "all",
                       !is.na(ID) & padj < 0.5 ~ "GWAS",
                       !is.na(ID) & padj >= 0.5 ~ "GWAS-NS")) %>%
  arrange(g) 

ggplot(comb_sort, aes(x= con, y = cold)) + 
  geom_point(aes(x = con, y = cold, color = g), alpha = 0.5) +
  scale_color_manual(values = c("#999999", "#E69F00", "#ff0000", "#000000"), breaks = "all") +
  scale_x_log10() +
  scale_y_log10() + 
  theme(legend.title = element_blank())


  

         