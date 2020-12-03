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

# Download and save GTF from flybase ---
# ftp://ftp.flybase.net/releases/FB2020_03/dmel_r6.34/gtf/
gtf_df <- read_clean_gtf(here::here("dmel-all-r6.34.gtf"))
# Add upstream and downstream annotation to it
gtf_df <- gtf_up_down_stream_annot(gtf_df)


# Group and filter snps -----
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




# Annotate the genes and features associated with each snp
fst_ch_snp_top <- gene_assoc_snp(fst_ch_snp_top, gtf_df)
saveRDS(fst_ch_snp_top, here::here("por/fst_ch_snp_top_annot.rds"))

x <- readRDS(here::here("por/fst_ch_snp_top_annot.rds"))





