# ------------------------------------------------------------------------------
# Por's Introgression data - Stats on allele freqs bw replicates and VT
# December 07, 2020
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
  
  
  # Pooling allele counts across groups -----
  
  # Create another nested data.frame from the allele counts that pools 
  # all replicates together to be compared to the VT background
  df[[i]]$allele_pooled_df <- map(df[[i]]$allele_count_df, 
                                  pool_alleles_all)
  
  # Create another nested data.frame from the allele counts that pools 
  # all replicates of the specific reciprocal cross -- ditches the VT parental
  df[[i]]$allele_rec_cross_pooled_df <- map(df[[i]]$allele_count_df, 
                                            pool_alleles_rec_cross)
  
  
  # Run the statistical tests ----
  
  # Individual pairwise Fisher tests sum-of-logs combined p-value
  df[[i]]$fish_pval_pw_comb <- map_dbl(df[[i]]$allele_count_df, 
                                       run_fisher_pw_comb)
  
  # Pooled: hybrid v parental
  df[[i]]$fish_pval_hy_par <- map_dbl(df[[i]]$allele_pooled_df, 
                                      run_fisher_exact)
  
  # Pooled: reciprocal cross (e.g. SKF v VT8F)
  df[[i]]$fish_pval_rcross <- map_dbl(df[[i]]$allele_rec_cross_pooled_df, 
                                      run_fisher_exact)

  tictoc::toc()
}


# Bind all rows back together
df <- bind_rows(df)


# Adjust the p-values with Benjamini-Hochberg correction -----------------------


df$fish_pval_pw_comb_padj <- p.adjust(df$fish_pval_pw_comb, 
                                       method = "BH", 
                                       n = nrow(df))

df$fish_pval_hy_par_padj <- p.adjust(df$fish_pval_hy_par, 
                                       method = "BH", 
                                       n = nrow(df))

df$fish_pval_rcross_padj <- p.adjust(df$fish_pval_rcross, 
                                       method = "BH", 
                                       n = nrow(df))


# Save final output as a rds file -----
saveRDS(df, here::here("por/introgression_df_fisher_whole_genome.rds"))
