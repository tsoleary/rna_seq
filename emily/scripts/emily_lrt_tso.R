# ------------------------------------------------------------------------------
# Emily's transcriptomics data
# September 12, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)
require(DESeq2)
source(here::here("functions.R"))

# Import and round data from Salmon
counts <- round(read.table(here::here("emily/counts/Dm_countsMatrix.txt"),
                           header = TRUE))

# Load the metadata
metadata <- read_delim(here::here("emily/counts/emily_metadata.txt"), 
                                  delim = "\t",
                       col_types = cols(
                         sampleID = col_character(),
                         pop = col_factor(),
                         temp = col_factor(),
                         rep = col_factor(),
                         region = col_factor()))

# LRT: full_model = temp + region + temp:region
# dds <- DESeqDataSetFromMatrix(countData = counts, 
#                               colData = metadata, 
#                               design = ~  temp + region + temp:region)
# Emily has shown that there was no significant results for the interaction 
# so I am going to ignore that for the sake of this script ...

# design without the interaction term 
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = ~  temp + region)

# LRT of Temperature
ddslrttemp <- DESeq(dds, test = "LRT",  reduced = ~ region)
temp_32_vs_25 <- results(ddslrttemp, name = "temp_32_vs_25", 
                         alpha = 0.05,
                         tidy = TRUE)
temp_34_vs_25 <- results(ddslrttemp, name = "temp_34_vs_25", 
                         alpha = 0.05,
                         tidy = TRUE)
temp_36_vs_25 <- results(ddslrttemp, name = "temp_36_vs_25", 
                         alpha = 0.05,
                         tidy = TRUE)

# LRT of Region
ddslrtreg <- DESeq(dds, test = "LRT",  reduced = ~ temp)
region_trop_vs_temp <- results(ddslrtreg, 
                               name = "region_temperate_vs_tropical", 
                               alpha = 0.05,
                               tidy = TRUE)

saveRDS(ddslrtreg, here::here("emily/counts/ddslrtreg.rds"))
saveRDS(region_trop_vs_temp, 
        here::here("emily/counts/region_trop_vs_temp.rds"))

x <- region_trop_vs_temp %>%
  filter(padj < 0.05)

gtf_df <- read_clean_gtf("~/Downloads/dmel-all-r6.34.gtf")
region_genes <- transcript_to_gene(x$row, gtf_df)

which(region_genes == "Ucp4A")
