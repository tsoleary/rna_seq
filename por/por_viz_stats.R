# ------------------------------------------------------------------------------
# Por's Introgression data - Visualization
# December 08, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages and functions
require(tidyverse)
require(qqman)

# Load data
df <- readRDS("por/introgression_df_fisher_whole_genome.rds")

chromosomes <- c("2L", "2R", "3L", "3R", "4", "X", "Y")

# Use a subset just to make trial plots faster
z <- df[sample(1:nrow(df), 100000, replace = FALSE), ] %>%
  filter(chr %in% chromosomes)

# Correlations ----

# Pairwise Fisher Exact Tests combined versus Pooled Replicate test
ggplot(z, aes(x = -log(fish_pval_pw_comb_padj), 
              y = -log(fish_pval_hy_par_padj))) +
  geom_point() +
  geom_smooth(method = "lm", 
              linetype = "dashed", 
              se = FALSE, formula = y ~ x, color = "grey50") +
  annotate("text", x = 60, y = 30, 
           label = paste("r =", 
                         round(cor(-log(z$fish_pval_pw_comb_padj), 
                                   -log(z$fish_pval_hy_par_padj)), 2)), 
           size = 4) + 
  theme_classic()

# Pooled Replicate Fisher Exact test versus Reciprocal cross
ggplot(z, aes(x = -log(fish_pval_hy_par_padj), 
              y = -log(fish_pval_rcross_padj))) +
  geom_point(aes(color = chr)) +
  geom_smooth(method = "lm", 
              linetype = "dashed", 
              se = FALSE, formula = y ~ x, color = "grey50") +
  annotate("text", x = 18, y = 12.5, 
           label = paste("r =", 
                         round(cor(-log(z$fish_pval_pw_comb_padj), 
                                   -log(z$fish_pval_rcross_padj)), 2)), 
           size = 4) + 
  theme_classic()


# Manhattan Plot 
z_manhat <- df %>%
  filter(chr %in% chromosomes) %>%
  mutate(CHR = case_when(chr == "2L" ~ 1,
                         chr == "2R" ~ 2,
                         chr == "3L" ~ 3,
                         chr == "3R" ~ 4,
                         chr == "4" ~ 5,
                         chr == "X" ~ 6,
                         chr == "Y" ~ 7))

manhattan(z_manhat, 
          chr = "CHR", 
          bp = "pos", 
          p = "arc_freq_ttest_padj", 
          chrlabs = chromosomes, 
          snp = "rc",
          suggestiveline = FALSE,
          genomewideline = -log10(1e-20))


