# Loading Emily's embryo transcriptomics data ----------------------------------

# Packages
require(tidyverse)
require(DESeq2)


# Load counts data -------------------------------------------------------------

# # Set working directory
# setwd(here::here("emily/counts"))
# 
# # Import and round data from Salmon
# counts <- round(read.table("Dm_countsMatrix.txt", header = TRUE))
# 
# # Load the metadata
# metadata <- read_delim("emily_metadata.txt", delim = "\t",
#                        col_types = cols(
#                          sampleID = col_character(),
#                          pop = col_factor(),
#                          temp = col_factor(),
#                          rep = col_factor(),
#                          region = col_factor()))
#   
# 
# # Run DESeq2 -------------------------------------------------------------------
# metadata <- metadata %>%
#   mutate(group = as.factor(paste(region, temp, sep = "_")))
# 
# 
# dds <- DESeqDataSetFromMatrix(countData = counts,
#                               colData = metadata,
#                               design = ~ group)
# 
# # Set the group to be compared against
# dds$group <- relevel(dds$group, ref = "temperate_36")
# 
# # Filter out genes with few reads
# dds <- dds[rowSums(counts(dds)) > 110]
# 
# # Run DESeq
# dds <- DESeq(dds)
# 
# # List the results you've generated
# resultsNames(dds)
# 
# res <- results(dds, 
#                contrast = c("group", "tropical_36", "temperate_36"), 
#                alpha = 0.05)
# 
# # Save an object to a file
# setwd(here::here("emily/results"))
# saveRDS(res, file = "trop_temp_36_results.rds")

# Load DESeq2 results ----------------------------------------------------------
setwd(here::here("emily/results"))

# Restore the object from the rds file that you created
res_25 <- readRDS("trop_temp_25_results.rds")
res_32 <- readRDS("trop_temp_32_results.rds")
res_34 <- readRDS("trop_temp_34_results.rds")
res_36 <- readRDS("trop_temp_36_results.rds")

res_temp_32 <- readRDS("temp_32_25_results.rds")
res_temp_34 <- readRDS("temp_34_25_results.rds")
res_temp_36 <- readRDS("temp_36_25_results.rds")

res_trop_32 <- readRDS("trop_32_25_results.rds")
res_trop_34 <- readRDS("trop_34_25_results.rds")
res_trop_36 <- readRDS("trop_36_25_results.rds")
  


# GTF downloaded from this website
# ftp://ftp.flybase.net/releases/FB2020_03/dmel_r6.34/gtf/

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
# Function: transcript_to_gene
# Description: convert a vector of transcripts to their gene symbols
# Inputs: character vector of transcripts and gtf data frame
# Outputs: gene

require(tidyverse)

transcript_to_gene <- function (trans_vec, gtf_dat){
  gtf_dat <- gtf_dat %>%
    distinct(transcript_id, .keep_all = TRUE) %>%
    filter(transcript_id != "")
  dat <- tibble::enframe(trans_vec, name = NULL, value = "transcript_id")
  dat$gene <- dat$transcript_id
  for (i in 1:nrow(dat)){
    print(i)
    temp <- which(dat$transcript_id[i] == gtf_dat$transcript_id, TRUE)
    if (length(temp) == 0) {
      gene_replace <- NA
    } else {
      gene_replace <- gtf_dat$gene_symbol[temp]
    }
    dat$gene <- gsub(dat$transcript_id[i],
                     gene_replace,
                     dat$gene)
  }
  dat_wo_na <- dat$gene[which(!is.na(dat$gene))]
  return(dat_wo_na)
}
# End function -----------------------------------------------------------------

