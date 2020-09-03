# ------------------------------------------------------------------------------
# Por's Introgression Mapping Gene Positions
# September 02, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)
source(here::here("functions.R"))

# Load the Fst tables ---
fst_ch <- read_csv(here::here("por/fst_CHxVT10_500_tidy.csv")) %>%
  mutate(start = BP - 250,
         end = BP + 250)
fst_sk <- read_csv(here::here("por/fst_SKxVT8_500_tidy.csv")) %>%
  mutate(start = BP - 250,
         end = BP + 250)

# SNP level Fst table 
fst_ch_snp <- read_csv(here::here("por/fst_CHxVT10_tidy.csv"))
fst_sk_snp <- read_csv(here::here("por/fst_SKxVT8_tidy.csv"))


# Download and save GTF from flybase ---
# ftp://ftp.flybase.net/releases/FB2020_03/dmel_r6.34/gtf/
gtf_df <- read_clean_gtf("~/Downloads/dmel-all-r6.34.gtf")

# Annotate 500 bp windows with overlapping genes
fst_ch <- gene_assoc_window(fst_ch, gtf_df)
fst_sk <- gene_assoc_window(fst_sk, gtf_df)

write_csv(fst_ch, here::here("por/fst_CHxVT10_500_tidy_gene_annot.csv"))
write_csv(fst_sk, here::here("por/fst_SKxVT8_500_tidy_gene_annot.csv"))




# Annotate snps with associated genes
# Probably should filter the top snps before doing the gene annotation because 
# otherwise it would take far too long.
fst_ch_snp_top <- fst_ch_snp %>%
  filter(CHvCHF > quantile(CHvCHF, 0.95)) %>%
  filter(!is.na(gene_assoc))

fst_ch_snp <- gene_assoc_window(fst_ch_snp, x)
fst_sk_snp <- gene_assoc_window(fst_sk_snp, x)

write_csv(fst_ch_snp, here::here("por/fst_CHxVT10_tidy_gene_annot.csv"))
write_csv(fst_sk_snp, here::here("por/fst_SKxVT8_tidy_gene_annot.csv"))


# Pick outliers based on a specific comparison for sliding window Fst tables
fst_ch_top <- fst_ch %>%
  filter(CHvCHF > quantile(CHvCHF, 0.95)) %>%
  filter(!is.na(gene_assoc))

fst_sk_top <- fst_sk %>%
  filter(CHvCHF > quantile(CHvCHF, 0.95)) %>%
  filter(!is.na(gene_assoc))
  

# Clean up gene annotation column for downstream analysis, functional enrichment
num_genes <- max(stringr::str_count(fst_ch_top$gene_assoc, pattern = ";") + 1)

fst_ch_top <- fst_ch_top %>%
  separate(gene_assoc, 
           into = paste("gene", 1:num_genes, sep = "_"), 
           sep = ";") %>%
  pivot_longer(contains("gene_"), 
               names_to = "lab", 
               values_to = "gene",
               values_drop_na = TRUE) %>% 
  distinct(gene)

# Save the file as a tab delim .txt file for webgesault
write_delim(fst_ch_top, 
            here::here("por/fst_ch_genes_CHvCHF_top_5_percent.txt"), 
            delim = "/t",
            col_names = FALSE)

