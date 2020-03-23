# Loading Emily's embryo transcriptomics data ----------------------------------

# Packages
require(tidyverse)
require(DESeq2)


# Load counts data -------------------------------------------------------------

# Set working directory
setwd(here::here("emily/counts"))

# Import and round data from Salmon
counts <- round(read.table("Dm_countsMatrix.txt", header = TRUE))

# Load the metadata
metadata <- read_delim("emily_metadata.txt", delim = "\t",
                       col_types = cols(
                         sampleID = col_character(),
                         pop = col_factor(),
                         temp = col_factor(),
                         rep = col_factor(),
                         region = col_factor()))
  

# Run DESeq2 -------------------------------------------------------------------
metadata <- metadata %>%
  mutate(group = as.factor(paste(region, temp, sep = "_")))


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ group)

# Set the group to be compared against
dds$group <- relevel(dds$group, ref = "tropical_25")

# Filter out genes with few reads
dds <- dds[rowSums(counts(dds)) > 110]

# Run DESeq
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

res <- results(dds, contrast = c("group", "tropical_25", "tropical_36"))

resdata <- merge(as.data.frame(res),
                 as.data.frame(counts(dds, normalized = TRUE)),
                 by = "row.names",
                 sort = FALSE)
names(resdata)[1] <- "gene"

# Write results to file
setwd(here::here("emily/results"))
write.csv(resdata,
          file = "trop_36_25_results_with_norm.csv",
          row.names = FALSE)

# Load DESeq2 results ----------------------------------------------------------
setwd(here::here("emily/results"))
res_25 <- read_csv("trop_temp_25_results_with_norm.csv")
res_32 <- read_csv("trop_temp_32_results_with_norm.csv")
res_34 <- read_csv("trop_temp_34_results_with_norm.csv")
res_36 <- read_csv("trop_temp_36_results_with_norm.csv")

res_temp_32 <- read_csv("temp_32_25_results_with_norm.csv")
res_temp_34 <- read_csv("temp_34_25_results_with_norm.csv")
res_temp_36 <- read_csv("temp_36_25_results_with_norm.csv")

res_trop_32 <- read_csv("trop_32_25_results_with_norm.csv")
res_trop_34 <- read_csv("trop_34_25_results_with_norm.csv")
res_trop_36 <- read_csv("trop_36_25_results_with_norm.csv")
  


# GTF ----- nonsense!
setwd("~/Downloads")

x <- read_delim("dmel-all-r6.32.gtf", 
                delim = "\t", 
                col_names = c("chr", "source", "feature", "start", 
                              "end", "score", "strand", "frame", "attributes"))

x$attributes <- str_replace_all(x$attributes, "\"", "")

x <- x %>%
  separate(attributes, 
           sep = ";", 
           into = c("gene_id", "gene_symbol", "transcript_id", 
                    "transcript_symbol", "notes", "extra", "extraextra"))
x$gene_id <- str_replace_all(x$gene_id, "gene_id ", "")
x$gene_symbol <- str_replace_all(x$gene_symbol, " gene_symbol ", "")
x$transcript_id <- str_replace_all(x$transcript_id, " transcript_id ", "")
x$transcript_symbol <- str_replace_all(x$transcript_symbol, " transcript_symbol ", "")


x <- x %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
  filter(transcript_id != "")


sum(res_25$gene %in% x$transcript_id)


# messing with this not done --

transcript_to_gene <- function (dat, gtf_dat){
  gtf_dat <- gtf_dat %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
    filter(transcript_id != "")
  
  dat$transcript_symbol <- dat$gene
  
  for (i in 1:nrow(dat)){
    temp <- which(dat$transcript_symbol[i] == gtf_dat$transcript_id, TRUE)
    dat$gene <- gsub(dat$transcript_symbol[i], 
                     gtf_dat$transcript_id[temp], 
                     dat$transcript_symbol)
  }
  
  return(dat$gene)
}
