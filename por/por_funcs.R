# ------------------------------------------------------------------------------
# Por's data quick functions
# December 02, 2020
# TS O'Leary
# ------------------------------------------------------------------------------



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
# Function: run_fisher
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