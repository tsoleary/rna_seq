# ------------------------------------------------------------------------------
# Oxidative stress gene results
# March 20, 2021
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)

ddslrtreg <- readRDS(here::here("emily/counts/ddslrtreg.rds"))
region_res <- read_csv(here::here("emily/results/tropvtempREGION.csv"))
genes_df <- read_csv(here::here("emily/results/emily_transcript_genes.csv"))

ox_stress_genes <- c("Sod3", "Ucp4A", "Atg6", "sesB", "rl", 
                     "Ccs", "Gyc89Db",  "Keap1", "CG6762", 
                     "Sdhaf3", "dj-1beta", "Prx2540-1", "Idh")

ox_stress_trans <- genes_df %>%
  filter(gene %in% ox_stress_genes)

x <- region_res %>%
  filter(transcript_id %in% ox_stress_trans$transcript_id,
         padj < 0.05, 
         log2FoldChange > 0) %>%
  left_join(ox_stress_trans) %>%
  arrange(padj) %>%
  group_by(gene) %>%
  #top_n(1, wt = baseMean) %>%
  select(transcript_id, gene, baseMean, log2FoldChange, padj) %>%
  mutate(padj = round(padj, digits = 5),
         baseMean = round(baseMean, digits = 2),
         log2FoldChange = round(log2FoldChange, digits = 2))


# Jill Question ----------------------------------------------------------------


# Data wrangling ----- 

# Load data from Fronteirs
df_f <- read_csv(here::here("cahan/results/Dm_cahan_deg_hot_results_with_norm.csv"))

df_32 <- read_csv(here::here("emily/results/temp_32_25_results_with_norm.csv")) %>%
  rename(transcript_id = gene) %>%
  left_join(genes_df)
df_34 <- read_csv(here::here("emily/results/temp_34_25_results_with_norm.csv")) %>%
  rename(transcript_id = gene) %>%
  left_join(genes_df)
df_36 <- read_csv(here::here("emily/results/temp_36_25_results_with_norm.csv")) %>%
  rename(transcript_id = gene) %>%
  left_join(genes_df)

# Row bind the dataframe
df <- bind_rows(list("adult" = df_f, 
                     "32" = df_32, 
                     "34" = df_34, 
                     "36" = df_36), 
                     .id = "heat")

# Total number of genes or transcripts -----------------------------------------

# Count the number of genes or transcripts expressed!
df %>%
  group_by(heat) %>%
  filter(!is.na(padj)) %>%
  count()

# Count the number of unique genes expressed
df %>%
  group_by(heat, gene) %>%
  top_n(n = 1, wt = baseMean) %>%
  filter(!is.na(padj)) %>%
  count()

# Total number of significantly up or down DEG(/Transcripts) in each -----------

# Filter only significant genes & transcripts
df_sig <- df %>%
  filter(padj < 0.05) 

# Count transcripts significantly up and down in each
df_sig %>%
  group_by(heat, log2FoldChange < 0) %>%
  count()

# Count unique genes significantly up or down in each
df_sig %>%
  group_by(heat, log2FoldChange < 0, gene) %>%
  top_n(n = 1, wt = baseMean) %>%
  count()

# Overlap between temp responsive genes in each category -----------------------

# Up-regulated genes
df_sig %>%
  group_by(heat, gene) %>%
  top_n(n = 1, wt = baseMean) %>%
  filter(log2FoldChange > 0) %>%
  count() %>%
  spread(heat, n, fill = 0) %>%
  ungroup(gene) %>%
  select(-gene) %>% 
  as.matrix() %>% 
  crossprod()

# Down-regulated genes
df_sig %>%
  group_by(heat, gene) %>%
  top_n(n = 1, wt = baseMean) %>%
  filter(log2FoldChange < 0) %>%
  count() %>%
  spread(heat, n, fill = 0) %>%
  ungroup(gene) %>%
  select(-gene) %>% 
  as.matrix() %>% 
  crossprod()

# Shared up-regulated gene list ------------------------------------------------

# Which genes are shared?
z <- df_sig %>%
  group_by(heat, gene)  %>%
  top_n(n = 1, wt = baseMean) %>%
  filter(log2FoldChange > 0) %>%
  group_by(gene) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>%
  select(gene, heat) %>%
  group_by(gene) %>%
  summarise(groups = paste(sort(unique(heat)), collapse = ", "))





