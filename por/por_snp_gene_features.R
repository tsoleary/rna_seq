# ------------------------------------------------------------------------------
# Annotate Por's SNPs to gene and genomic features
# October 20, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)
source(here::here("functions.R"))

# SNP level Fst table 
fst_ch_snp <- read_csv(here::here("por/fst_CHxVT10_tidy.csv"))
#fst_sk_snp <- read_csv(here::here("por/fst_SKxVT8_tidy.csv"))


# Download and save GTF from flybase ---
# ftp://ftp.flybase.net/releases/FB2020_03/dmel_r6.34/gtf/
gtf_df <- read_clean_gtf(here::here("dmel-all-r6.34.gtf"))

# Annotate snps with associated genes 

# Group and filter snps

# Proportion of top snps to keep
p_top_snps <- 0.1
quant_p <- 1 - p_top_snps

fst_ch_snp_top <- fst_ch_snp %>%
  mutate(group = case_when(VT10vVT10F > quantile(VT10vVT10F, quant_p) &
           CHFvVT10 > quantile(CHFvVT10, quant_p) ~ "VT10vVT10F_CHFvVT10",
         VT10vVT10F > quantile(VT10vVT10F, quant_p) & 
           !(CHFvVT10 > quantile(CHFvVT10, quant_p)) ~ "VT10vVT10F",
         !(VT10vVT10F > quantile(VT10vVT10F, quant_p)) & 
           CHFvVT10 > quantile(CHFvVT10, quant_p) ~ "CHFvVT10")) %>%
  filter(!is.na(group))


# fst_sk_snp_top <- fst_sk_snp %>%
#   mutate(group = case_when(SKFvVT8 > quantile(SKFvVT8, quant_p) &
#            VT8vVT8F > quantile(VT8vVT8F, quant_p) ~ "SKFvVT8_VT8vVT8F",
#          SKFvVT8 > quantile(SKFvVT8, quant_p) &
#            !(VT8vVT8F > quantile(VT8vVT8F, quant_p)) ~ "SKFvVT8",
#          !(SKFvVT8 > quantile(SKFvVT8, quant_p)) &
#            VT8vVT8F > quantile(VT8vVT8F, quant_p) ~ "VT8vVT8F")) %>%
#   filter(!is.na(group))

gtf_df <- gtf_up_down_stream_annot(gtf_df)

fst_ch_snp_top <- gene_assoc_snp_pmap(fst_ch_snp_top, gtf_df)
#fst_sk_snp_top <- gene_assoc_snp(fst_sk_snp_top, gtf_df)




write_csv(fst_ch_snp, here::here("por/fst_CHxVT10_tidy_gene_annot.csv"))
#write_csv(fst_sk_snp, here::here("por/fst_SKxVT8_tidy_gene_annot.csv"))


# Clean up gene annotation column for downstream analysis, functional enrichment
num_genes <- max(stringr::str_count(fst_ch_top$gene_assoc, pattern = "|") + 1)

fst_ch_top <- fst_ch_top %>%
  separate(gene_assoc, 
           into = paste("gene", 1:num_genes, sep = "_"), 
           sep = "|") %>%
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





