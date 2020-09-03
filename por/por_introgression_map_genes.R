# ------------------------------------------------------------------------------
# Por's Introgression Mapping Gene Positions
# September 02, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)

# Load the Fst tables ---
fst_ch <- read_csv(here::here("por/fst_CHxVT10_500_tidy.csv")) %>%
  mutate(start = BP - 250,
         end = BP + 250)
fst_sk <- read_csv(here::here("por/fst_SKxVT8_500_tidy.csv")) %>%
  mutate(start = BP - 250,
         end = BP + 250)


# Download and save GTF from flybase ---
# ftp://ftp.flybase.net/releases/FB2020_03/dmel_r6.34/gtf/
x <- read_clean_gtf("~/Downloads/dmel-all-r6.34.gtf")


# Function to match them up


# ------------------------------------------------------------------------------
# Function: gene_assoc_window
# Description: Annotate with all genes associated with the specific window
# Inputs: fst data frame and gtf data frame
# Outputs: fst data frame with new gene_assoc column

require(tidyverse)

gene_assoc_window <- function(dat, gtf_dat) {
  
  dat$gene_assoc <- vector(mode = "character", length = nrow(dat))
  
  gtf_dat <- gtf_dat %>%
    filter(feature == "gene")
  
  
  for (i in 1:nrow(dat)){
    print(i)
    
    # Find all genes that overlap the window in anyway
    temp <- c( 
      which(dat$CHR[i] == gtf_dat$chr & 
              gtf_dat$end > dat$start[i] & 
              gtf_dat$end < dat$end[i], TRUE),
      which(dat$CHR[i] == gtf_dat$chr & 
              gtf_dat$start > dat$start[i] & 
              gtf_dat$start < dat$end[i], TRUE),
      which(dat$CHR[i] == gtf_dat$chr & 
              gtf_dat$start < dat$start[i] & 
              gtf_dat$end > dat$end[i], TRUE)
    )
    
    # Paste all gene symbols together with ; between. NA if there is none
    if (is_empty(temp)){
      dat$gene_assoc[i] <- NA
    } else{
      dat$gene_assoc[i] <- paste(unique(gtf_dat$gene_symbol[temp]), 
                                 collapse = ";")
    }
    
  }
  return(dat)
} 
# End function -----------------------------------------------------------------


fst_ch <- gene_assoc_window(fst_ch, x)
fst_sk <- gene_assoc_window(fst_sk, x)





