# ------------------------------------------------------------------------------
# Por's Introgression data - Stats on allele freqs bw replicates and VT
# December 02, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages and functions
require(tidyverse)
source(here::here("por/por_funcs.R"))

# Load data
df_messy <- read_delim("~/Downloads/SKxVT8_3rep_rc", 
                       delim = "\t") %>% head(1000)
# Running this on the whole thing currently just tanks it. 
# So I think either for looping each row individually or like a thousand 
# 1000 row chunks make work


# Tidy the data
df <- df_messy %>%
  # Pivot the data frame to long format to separate out type of allele 
  # and count from total
  pivot_longer(contains(c("_maj", "_min")), 
               names_to = "geno_allele", 
               values_to = "allele_frac") %>%
  # Separate out the "genotype" from the allele typle
  separate(geno_allele, 
           into = c("geno", "major_minor"), 
           sep = "_") %>%
  # Remove SK parent bc we won't use it down stream
  filter(geno != "SK") %>%
  # Separate out the specific allele count from the total
  separate(allele_frac, 
           into = c("count", "total"), 
           sep = "/") %>%
  # Remove the total because we don't need it
  select(-total) %>%
  # Need to treat the counts as a numeric 
  # (it is currently a character class bc the original "97/100" format)
  mutate(count = as.numeric(count)) %>%
  # Flip it back wider so that genotype, major, and minor allele counts have
  # their own column each
  pivot_wider(values_from = "count", 
              names_from = "major_minor") %>%
  # Nest all that in a data.frame so that each row is back to a single loci!
  nest(allele_count_df = c(geno, maj, min))

# Pooling allele counts across groups ------------------------------------------

# Create another nested data.frame from the allele counts that pools 
# all replicates together to be compared to the VT background
df$allele_pooled_df <- map(df$allele_count_df, 
                           pool_alleles_all)

# Create another nested data.frame from the allele counts that pools 
# all replicates of the specific reciprocal cross -- ditches the VT parental
df$allele_rec_cross_pooled_df <- map(df$allele_count_df, 
                                     pool_alleles_rec_cross)


# Run the two-sided Fisher exact test on each of the different comparisons -----

# Each individually
df$fish_pval_all <- map_dbl(df$allele_count_df, 
                            run_fisher)

df$fish_pval_all <- map_dbl(df$allele_count_df, 
                            run_fisher_MC)


for (i in 1:40){
  print(i)
  tictoc::tic()
  df$fisher_for_exact[i] <- run_fisher_exact(df$allele_count_df[[i]])
  tictoc::toc()
  tictoc::tic()
  df$fisher_for_MC_10000[i] <- run_fisher_MC(df$allele_count_df[[i]], 
                                         reps = 10000)
  tictoc::toc()
  tictoc::tic()
  df$fisher_for_MC_100000[i] <- run_fisher_MC(df$allele_count_df[[i]], 
                                         reps = 100000)
  tictoc::toc()
  tictoc::tic()
  df$fisher_for_MC_1000000[i] <- run_fisher_MC(df$allele_count_df[[i]], 
                                         reps = 1000000)
  tictoc::toc()
  tictoc::tic()
  df$fisher_for_MC_10000000[i] <- run_fisher_MC(df$allele_count_df[[i]], 
                                               reps = 10000000)
  tictoc::toc()
}

x <- df[1:40, ] %>% 
  select(chr, pos, contains("for_exact"), 21:24) %>%
  pivot_longer(contains("for_"), 
               names_to = "test", 
               values_to = "pval") %>%
  group_by(chr, pos) %>%
  mutate(delta_pval = pval - pval[test == "fisher_for_exact"]) %>%
  ungroup(chr, pos) %>%
  group_by(test) %>%
  summarize(sum_abs_delta_pval = sum(abs(delta_pval), na.rm = TRUE)) %>%
  mutate(pval_ind_off = sum_abs_delta_pval/40)



# Pooled: hybrid v parental
df$fish_pval_hy_par <- map_dbl(df$allele_pooled_df, 
                               run_fisher)

# Pooled: reciprocal cross (e.g. SKF v VT8F)
df$fish_pval_rcross <- map_dbl(df$allele_rec_cross_pooled_df, 
                               run_fisher)

write_delim(df %>% 
              select(-contains("df")), 
            here::here("por/results_temp_1000.txt"))


