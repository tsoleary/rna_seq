# ------------------------------------------------------------------------------
# Por's Introgression data - Stats on allele freqs bw replicates and VT
# January 25, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages and functions
require(tidyverse)
source(here::here("por/por_funcs.R"))

# Load data
df_messy <- read_delim("~/Downloads/SKxVT8_3rep_rc", 
                       delim = "\t")

# Split it into chunks
df <- chunkify(df_messy, r_chunk = 1000)

# Tidy the data and create other pooled allele count nested dfs and run tests --
for (i in 1:length(df)){
  print(i)
  tictoc::tic()
  
  # Tidy the data -----
  df[[i]] <- tidy_allele_df(df[[i]])

  # Calculate the p-value for all the One sample t-tests on arcsine transformed 
  # allele frequencies
  df[[i]]$arc_freq_ttest_pval <- map_dbl(df[[i]]$allele_count_df, 
                                         arc_freq_ttest)
  
  tictoc::toc()
}


# Bind all rows back together
df <- bind_rows(df)


# Adjust the p-values with Benjamini-Hochberg correction -----------------------

df <- df %>%
        filter(!is.na(arc_freq_ttest_pval))

df$arc_freq_ttest_padj <- p.adjust(df$arc_freq_ttest_pval, 
                                   method = "BH", 
                                   n = nrow(df))



# Save final output as a rds file -----
saveRDS(df, here::here("por/introgression_SK_VT8_25JAN2021_na_rm.rds"))
