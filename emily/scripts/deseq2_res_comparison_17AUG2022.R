# ------------------------------------------------------------------------------
# Comparing DESeq2 results with LRT and Wald
# August 17, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)
require(DESeq2)
#source(here::here("functions.R"))

# Import and round data from Salmon
counts <- round(read.table(here::here("emily/counts/Dm_countsMatrix.txt"),
                           header = TRUE))

# Load the metadata
metadata <- read_delim(here::here("emily/counts/metadata.txt"), 
                                  delim = "\t",
                       col_types = cols(
                         sampleID = col_character(),
                         pop = col_factor(),
                         temp = col_factor(),
                         rep = col_factor(),
                         region = col_factor()))

# LRT: full_model = temp + region
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ temp + region)

# LRT of Temperature
dds_lrt_temp <- DESeq(dds, test = "LRT",  reduced = ~ region)
res_32 <- results(dds_lrt_temp, name = "temp_32_vs_25", 
                  alpha = 0.05,
                  tidy = TRUE)
res_34 <- results(dds_lrt_temp, name = "temp_34_vs_25", 
                  alpha = 0.05,
                  tidy = TRUE)
res_36 <- results(dds_lrt_temp, name = "temp_36_vs_25", 
                  alpha = 0.05,
                  tidy = TRUE)

res_lrt_temp <- bind_rows(list("32" = res_32,
                               "34" = res_34,
                               "36" = res_36),
                          .id = "hs_temp")

# Wald of Temperature
dds_wald_temp <- DESeq(dds, test = "Wald")
res_32 <- results(dds_wald_temp, name = "temp_32_vs_25", 
                  alpha = 0.05,
                  tidy = TRUE)
res_34 <- results(dds_wald_temp, name = "temp_34_vs_25", 
                  alpha = 0.05,
                  tidy = TRUE)
res_36 <- results(dds_wald_temp, name = "temp_36_vs_25", 
                  alpha = 0.05,
                  tidy = TRUE)

res_wald_temp <- bind_rows(list("32" = res_32,
                                "34" = res_34,
                                "36" = res_36),
                           .id = "hs_temp")

res_temp <- bind_rows(list("Wald" = res_wald_temp,
                           "LRT" = res_lrt_temp),
                      .id = "Test")


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
