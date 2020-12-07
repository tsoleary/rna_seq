# ------------------------------------------------------------------------------
# Por's data quick functions
# December 02, 2020
# TS O'Leary
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: chunkify
# Description: split a data.frame to a list of data.frames of workable sizes
# Inputs: the big data.frame and r_chunk is the rows per chunk
# Outputs: list of data.frames

chunkify <- function(df, r_chunk = 1000) {
  
  chunk  <- rep(1:ceiling(nrow(df)/r_chunk), each = r_chunk)[1:nrow(df)]
  
  return(split(df, chunk))
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: tidy_allele_df
# Description: tidy the allele count df into nested df for fisher tests
# Inputs: data.frame from Por's allele count output
# Outputs: tidy nested data.frame for fisher tests

tidy_allele_df <- function(dat) {
  dat %>%
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
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: pool_alleles_all
# Description: Pool the all replicate hybrid alleles compared to VT background
# Inputs: data.frame with a geno, min, and maj columns
# Outputs: data frame with allele counts pooled

require(tidyverse)

pool_alleles_all <- function(dat) {
  dat %>% 
    mutate(geno = case_when(str_detect(geno, "F") == TRUE ~ "hybrid",
                            str_detect(geno, "F") == FALSE ~ "parental")) %>%
    group_by(geno) %>%
    summarize(maj = sum(maj), 
              min = sum(min), 
              .groups = "drop")
}
# End function -----------------------------------------------------------------



# ------------------------------------------------------------------------------
# Function: pool_alleles_rec_cross
# Description: Pool the replicate reciprocal crosses to compare each other
# Inputs: data.frame with a geno, min, and maj columns
# Outputs: data frame with allele counts pooled

require(tidyverse)

pool_alleles_rec_cross <- function(dat) {
  dat %>% 
    mutate(geno = case_when(str_detect(geno, "SKF") == TRUE ~ "SKF",
                            str_detect(geno, "VT8F") == TRUE ~ "VT8F")) %>%
    filter(!is.na(geno)) %>%
    group_by(geno) %>%
    summarize(maj = sum(maj), 
              min = sum(min), 
              .groups = "drop")
}
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: df_to_mat
# Description: quickly convert df to mat and set geno col to rownames
# Inputs: data.frame with a geno column
# Outputs: matrix with geno as the rownames

df_to_mat <- function(dat) {
  m <- as.matrix(dat[, -1])
  rownames(m) <- dat$geno
  return(m)
}
# End function -----------------------------------------------------------------




# ------------------------------------------------------------------------------
# Function: run_fisher_exact
# Description: run fisher exact test on a data.frame 
# Inputs: data.frame with a geno column
# Outputs: matrix with geno as the rownames

run_fisher_exact <- function(dat) {
  p_val <- tryCatch(fisher.test(df_to_mat(dat), 
                                alternative = "two.sided", 
                                workspace = 2000000)$p.value, 
                    error = function(err) NA)
  
  return(p_val)
}
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: run_fisher_MC
# Description: run fisher exact test on a data.frame 
# Inputs: data.frame with a geno column
# Outputs: matrix with geno as the rownames

run_fisher_MC <- function(dat, reps = 1000000) {
  p_val <- tryCatch(fisher.test(df_to_mat(dat), 
                                alternative = "two.sided", 
                                simulate.p.value = TRUE,
                                B = reps)$p.value, 
                    error = function(err) NA)
  
  return(p_val)
}
# End function -----------------------------------------------------------------



# ------------------------------------------------------------------------------
# Function: run_fisher_pw_comb
# Description: Run pairwise Fisher tests and combine BH corrected pvals with 
#              the sum of logs method
# Inputs: data.frame with geno column with the maj and min alleles, 
#         and a background to compare it to
# Outputs: single sum-of-logs combined p-value

run_fisher_pw_comb <- function(dat, background = "VT8") {
  
  # Get the rows of the data frame that match and don't match the background 
  # allele count that you are comparing each to
  r_background <- which(dat$geno == background)
  r_other <- which(dat$geno != background)
  
  # Initialize a vector to store the pairwise Fisher test p-values
  pvals <- vector(mode = "numeric", length = length(r_other))
  
  # Run the Fisher test on all pairwise 
  for (i in 1:length(r_other)) {
    pvals[i] <- run_fisher_exact(dat[c(r_background, r_other[i]), ])
  }
  
  # Adjust the p-values with Benjamini-Hochberg correction
  fdr <- p.adjust(pvals, method = "BH", n = length(pvals))
  
  # Combine p-values with the sum of logs method
  p_comb <- metap::sumlog(fdr)$p
  
  return(p_comb)
} 
# End function -----------------------------------------------------------------
