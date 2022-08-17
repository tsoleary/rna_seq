# ------------------------------------------------------------------------------
# Functions for gene annotations
# March 18, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Function: read_clean_gtf
# Description: Read in and clean gtf file
# Inputs: gtf_file
# Outputs: gtf cleaned data frame

require(tidyverse)

read_clean_gtf <- function(file, bp_flanking = 1000) {
  # Read in the gtf
  x <- read_delim(file, 
                  delim = "\t", 
                  col_names = c("chr", "source", "feature", "start", "end", 
                                "score", "strand", "frame", "attributes"))
  
  # Clean up and separate the attributes file
  x$attributes <- str_replace_all(x$attributes, "\"", "")
  x <- x %>%
    separate(attributes, 
             sep = ";", 
             into = c("gene_id", "gene_symbol", "transcript_id", 
                      "transcript_symbol", "notes", "extra", "extraextra"))
  # remove the name and the awkward space before
  x$gene_id <- str_replace_all(x$gene_id, 
                               "gene_id ", "")
  x$gene_symbol <- str_replace_all(x$gene_symbol, 
                                   " gene_symbol ", "")
  x$transcript_id <- str_replace_all(x$transcript_id, 
                                     " transcript_id ", "")
  x$transcript_symbol <- str_replace_all(x$transcript_symbol, 
                                         " transcript_symbol ", "")
  
  dat <- gtf_up_down_stream_annot(x, bp_up = bp_flanking)
  return(dat)
} 
# End function -----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Function: gene_assoc_snp
# Description: Annotate with all genes associated with the snp
# Inputs: data frame and gtf data frame
# Outputs: data frame with new gene_assoc column

require(tidyverse)

gene_assoc_snp <- function(dat, gtf_dat) {
  
  # Initialize a column to store the gene associations
  dat$gene_assoc <- vector(mode = "character", length = nrow(dat))
  
  # Loop through each row of the dataframe individually
  pb <- progress::progress_bar$new(total = nrow(dat))
  
  for (i in 1:nrow(dat)){
    # Print the row just to watch the progress
    pb$tick()
    
    # Find all genes that overlap the window in anyway
    temp <- c( 
      which(dat$CHROM[i] == gtf_dat$chr &  
              gtf_dat$start <= dat$POS[i] & 
              gtf_dat$end >= dat$POS[i], TRUE)
    )
    
    # Paste all gene symbols together with ; between. NA if there is none
    if (is_empty(temp)){
      dat$gene_assoc[i] <- NA
    } else{
      x <- gtf_dat[temp, ]
      
      # Split the df based on gene_symbol
      y <- x %>% 
        group_split(gene_symbol)
      
      # Initialize an empty object to store the string
      genes_feature_str <- NULL
      
      # Loop through each gene
      for (j in length(y)){
        # Get the gene and feature information. Features separated by a comma
        gene_str <- unique(y[[j]]$gene_symbol)
        # maybe add something to get the heirachy of features later
        feature_str <- paste(unique(y[[j]]$feature), 
                             collapse = ",") 
        # Paste the gene symbol to the associated features sep by semicolon
        gene_feature_str <- paste(gene_str, feature_str, sep = ";")
        
        genes_feature_str <- c(genes_feature_str, gene_feature_str)
      }
      
      dat$gene_assoc[i] <- paste(genes_feature_str, 
                                 collapse = "|")
      
    }
    
  }
  return(dat)
  
} 
# End function -----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Function: gtf_up_down_stream_annot
# Description: Adds upstream and downstream flanking sequence
# Inputs: gtf file
# Outputs: output_description

require(tidyverse)

gtf_up_down_stream_annot <- function(gtf_df, bp_up = 1000, bp_down = bp_up) {
  
  gtf_up_pos <- gtf_df %>%
    filter(feature == "gene", strand == "+") %>%
    mutate(end = start - 1,
           start = start - bp_up - 1,
           feature = "upstream")
  
  gtf_down_pos <- gtf_df %>%
    filter(feature == "gene", strand == "+") %>%
    mutate(start = end + 1,
           end = end + bp_down + 1,
           feature = "downstream")
  
  gtf_up_neg <- gtf_df %>%
    filter(feature == "gene", strand == "-") %>%
    mutate(end = start - 1,
           start = start - bp_down - 1,
           feature = "downstream")
  
  gtf_down_neg <- gtf_df %>%
    filter(feature == "gene", strand == "-") %>%
    mutate(start = end + 1,
           end = end + bp_up + 1,
           feature = "upstream")
  
  gtf_df <- bind_rows(gtf_down_pos, 
                      gtf_up_pos,
                      gtf_down_neg, 
                      gtf_up_neg,
                      gtf_df) %>% 
    arrange(desc(chr), start)
  
  return(gtf_df)
} 

# End function -----------------------------------------------------------------