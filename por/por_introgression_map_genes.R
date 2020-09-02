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
x <- read_delim("~/Downloads/dmel-all-r6.34.gtf", 
                delim = "\t", 
                col_names = c("chr", "source", "feature", "start", 
                              "end", "score", "strand", "frame", "attributes"))

x$attributes <- str_replace_all(x$attributes, "\"", "")

x <- x %>%
  separate(attributes, 
           sep = ";", 
           into = c("gene_id", "gene_symbol", "transcript_id", 
                    "transcript_symbol", "notes", "extra", "extraextra"))
x$gene_id <- str_replace_all(x$gene_id,
                             "gene_id ", "")
x$gene_symbol <- str_replace_all(x$gene_symbol, 
                                 " gene_symbol ", "")
x$transcript_id <- str_replace_all(x$transcript_id, 
                                   " transcript_id ", "")
x$transcript_symbol <- str_replace_all(x$transcript_symbol, 
                                       " transcript_symbol ", "")


# Try to match them up
sum(fst_ch$start > x$start & fst_ch$CHR == x$chr)

dat <- fst_ch

dat$gene_assoc <- vector(mode = "character", length = nrow(dat))
nrow(dat)
for (i in 1:nrow(dat)){
  chr <- dat$CHR[i]
  start <- dat$start[i]
  end <- dat$end[i]
  
  temp <- which(chr == x$chr & start > x$start & end < x$end, TRUE)
  
  if (is_empty(temp)){
    dat$gene_assoc[i] <- NA
  } else{
    dat$gene_assoc[i] <- paste(unique(x$gene_symbol[temp]), collapse = ";")
  }
  


}


