# ------------------------------------------------------------------------------
# Por's Introgression data - Stats on allele freqs bw replicates and VT
# January 26, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages and functions
require(tidyverse)
source(here::here("por/por_funcs.R"))

# Load data
df_messy <- read_delim("~/Downloads/VT_allele_count_SKxVT8.csv",
                       delim = ",")

colnames(df_messy) <- c("X1", "CHR", "BP", "sk.maj.allele", "vt8.maj.allele", 
                        "sk_sk.maj.count", "vt8_vt.count", "sk_cov", "vt8_cov",
                        "skf1_cov", "skf2_cov", "skf3_cov", "vt8f1_cov", 
                        "vt8f2_cov", "vt8f3_cov", "skf1_vt.count",
                        "skf2_vt.count", "skf3_vt.count", "vt8f1_vt.count",
                        "vt8f2_vt.count", "vt8f3_vt.count")

# Split it into chunks
df <- chunkify(df_messy, r_chunk = 1000)

# Tidy the data and create other pooled allele count nested dfs and run tests --
for (i in 1:length(df)){
  print(i)
  tictoc::tic()
  # Tidy the data and calculate the VT major allele frequency
  df[[i]] <- tidy_vt_alleles(df[[i]])
  tictoc::toc()
}

# Bind all rows back together
df <- bind_rows(df)


# Get a quick histogram of the parental lines VT major allele frequency
x <- df %>%
  select(CHR, BP, sk.maj.allele, vt8.maj.allele, geno, cov, vt_freq) %>%
  pivot_wider(names_from = "geno", values_from = c("vt_freq", "cov")) %>%
  filter(vt_freq_sk == 1 & vt_freq_vt8 == 1) 

x %>%
  ggplot() +
  geom_density(aes(x = vt_freq, 
                     fill = geno), 
                 bins = 10, 
                 color = "grey50", 
                 position = "dodge") + 
  theme_classic()


# Stats ------------------------------------------------------------------------


# Nest the data for analysis
df <- df %>% 
  nest(nested_df = c(geno, cov, vt.count, vt_freq))

# Split it into chunks
df <- chunkify(df, r_chunk = 1000)

# Loop through all chunks and run the One-Sample t-test
for (i in 1:length(df)){
  print(i)
  tictoc::tic()
  # Calculate the p-value for all the One sample t-tests on arcsine transformed 
  # allele frequencies
  df[[i]]$arc_freq_ttest_pval <- map_dbl(df[[i]]$nested_df, 
                                         arc_freq_ttest)
  tictoc::toc()
}


# Bind all rows back together
df <- bind_rows(df)


# Adjust p-values using Benjamini-Hochberg FDR correction
df$arc_freq_ttest_padj_BH <- p.adjust(df$arc_freq_ttest_pval, 
                                      method = "BH", 
                                      n = nrow(df))

# Adjust p-values using Bonferroni's correction
df$arc_freq_ttest_padj_bonf <- p.adjust(df$arc_freq_ttest_pval, 
                                        method = "bonferroni", 
                                        n = nrow(df))


# Save in wider form 
df_wider <- df %>% 
  unnest(nested_df) %>%
  pivot_wider(names_from = geno, 
              values_from = c(cov, vt.count, vt_freq))

df_wider_avg <- df %>% 
  unnest(nested_df) %>% 
  mutate(type = case_when(str_detect(geno, "f") == TRUE ~ "introgression",
                          str_detect(geno, "f") == FALSE ~ "parental")) %>%
  group_by(CHR, BP, type) %>%
  summarize(vt_freq_avg = mean(vt_freq),
            cov_avg = mean(cov)) %>%
  pivot_wider(names_from = type, 
              values_from = c(vt_freq_avg, cov_avg))

df_wider <- cbind(df_wider, df_wider_avg[, 3:6])


df_bonf_05 <- df_wider %>%
  filter(arc_freq_ttest_padj_bonf < 0.05) 

df_na <- df_wider %>% 
  filter(is.na(arc_freq_ttest_pval))

####
df_wider %>%
  head(100000) %>%
  ggplot(aes(x = vt_freq_avg_introgression, 
             y = vt_freq_vt8,
             color = arc_freq_ttest_padj_bonf < 0.05)) +
  geom_point() + 
  theme_classic()


df_bonf_05 %>%
  ggplot(aes(x = vt_freq_avg_introgression, 
             y = vt_freq_vt8)) +
  geom_point() + 
  theme_classic()


## Annotate
source(here::here("functions.R"))
gtf_df <- read_clean_gtf(here::here("dmel-all-r6.34.gtf"))


x <- gene_assoc_snp(df_bonf_05, gtf_df)       


x %>%
  select(gene_assoc, arc_freq_ttest_padj_bonf, vt_freq_vt8, vt_freq_avg_introgression) %>%
  arrange(vt_freq_avg_introgression)

y <- df_wider %>%
  select(arc_freq_ttest_padj_bonf, vt_freq_vt8, vt_freq_avg_introgression) %>%
  arrange(vt_freq_avg_introgression) %>%
  head(10000)

