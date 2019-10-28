# read in .cnt files and compile a table to be imported into PIVOT -------------
library(dplyr)

# sample, group, filename information dataframe
samples <- read.delim("/Users/tsoleary/R/rna_seq/DM6_counts/target.txt")

# working directory containg the .cnt files
setwd("/Users/tsoleary/R/rna_seq/DM6_counts")

# read in each .cnt and compile a dataframe of read counts for each sample
x <- data.frame()

for (i in 1:length(as.character(samples$files))){
  t <- read.table(as.character(samples$files)[i])
  if (i == 1){
    x <- t 
    colnames(x) <- c("gene", as.character(samples$label)[i])
  } else{
    colnames(t) <- c("gene", as.character(samples$label)[i])
    x <- full_join(x, t, by = "gene")
  }
}

# create .csv file to upload to PIVOT for DESeq normalization
write.csv(x, "dm_rna_seq_counts.csv", row.names = FALSE)

# load PIVOT -------------------------------------------------------------------
library(PIVOT)
pivot()
