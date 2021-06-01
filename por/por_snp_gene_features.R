# ------------------------------------------------------------------------------
# Annotate Por's SNPs to gene and genomic features
# June 1, 2021
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)

# Load functions -----

# ------------------------------------------------------------------------------
# Function: read_clean_gtf
# Description: Read in and clean gtf file
# Inputs: gtf_file
# Outputs: gtf cleaned data frame

read_clean_gtf <- function(file) {
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
  
  return(x)
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




# Load data
df_vt8 <- read_csv("por/vt8x_hampel_sigsnps.csv")
df_vt10 <- read_csv("por/vt10x_hampel_sigsnps.csv")
df_overlap <- read_csv("por/overlapSNP.csv")


# Load most current gtf genome annotation with corresponding release ---
# ftp://ftp.flybase.org/genomes/dmel/current/gtf/
gtf_df <- read_clean_gtf("dmel-all-r6.39.gtf")

# Annotate the genes and features associated with each snp
df_vt8 <- gene_assoc_snp(df_vt8, gtf_df)
df_vt10 <- gene_assoc_snp(df_vt10, gtf_df)
df_overlap <- gene_assoc_snp(df_overlap, gtf_df)

# Clean up the gene annotations
df_vt8 <- df_vt8 %>%
  separate(gene_assoc, into = c("gene", "info"), sep = ";") 

df_vt10 <- df_vt10 %>%
  separate(gene_assoc, into = c("gene", "info"), sep = ";") 

df_overlap <- df_overlap %>%
  separate(gene_assoc, into = c("gene", "info"), sep = ";") 

# Save annotations
write_csv(df_vt8, "por/vt8x_hampel_sigsnps_genes.csv")
write_csv(df_vt10, "por/vt10x_hampel_sigsnps_genes.csv")
write_csv(df_overlap, "por/overlapSNPS_genes.csv")

