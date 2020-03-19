# Loading Emily's embryo transcriptomics data ----------------------------------

# Packages
require(tidyverse)
require(DESeq2)


# Run DESeq2 -------------------------------------------------------------------

# Set wd
setwd(here::here("emily/counts"))

# Import and round data from Salmon
counts <- round(read.table("Dm_countsMatrix.txt", header = TRUE))


df <- data.frame(sampleID = colnames(counts))
df <- df %>%
  separate(sampleID, sep = "_", 
           into = c("pop", "temp", "rep"), remove = FALSE) %>%
  mutate(region = case_when(
    pop %in% c("BO", "CP", "GH", "GM", "SK") ~ "tropical",
    !(pop %in% c("BO", "CP", "GH", "GM", "SK")) ~ "temperate"))

write_tsv(df, "emily_metadata.txt")
  



