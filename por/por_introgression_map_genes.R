# ------------------------------------------------------------------------------
# Por's Introgression Mapping Gene Positions
# September 02, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)

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
x <- read_clean_gtf("~/Downloads/dmel-all-r6.34.gtf")


# Annotate 500 bp windows with overlapping genes
fst_ch <- gene_assoc_window(fst_ch, x)
fst_sk <- gene_assoc_window(fst_sk, x)

write_csv(fst_ch, here::here("por/fst_CHxVT10_500_tidy_gene_annot.csv"))
write_csv(fst_sk, here::here("por/fst_SKxVT8_500_tidy_gene_annot.csv"))



# Annotate snps with associated genes
fst_ch_snp <- gene_assoc_window(fst_ch_snp, x)
fst_sk_snp <- gene_assoc_window(fst_sk_snp, x)

write_csv(fst_ch_snp, here::here("por/fst_CHxVT10_tidy_gene_annot.csv"))
write_csv(fst_sk_snp, here::here("por/fst_SKxVT8_tidy_gene_annot.csv"))
