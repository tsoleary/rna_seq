# gwas dgrp --------------------------------------------------------------------

library(dplyr)
library(stringr)

max_gwas_top <- read.table(
  "/Users/tsoleary/R/rna_seq/DRGP_GWAS/CTmax/gwas.top.annot", header = TRUE)

min_gwas_top <- read.table(
  "/Users/tsoleary/R/rna_seq/DRGP_GWAS/CTmin/gwas.top.annot", header = TRUE)

snps <- min_gwas_top$ID

# fly base requires the SNPs in a certain format
chromosome <-  str_extract(snps, "[[:digit:]]?[[:alpha:]]")
position <- str_replace(snps, "[[:digit:]]?[[:alpha:]]_", "") %>%
  str_replace("_[[:alnum:]]+", "")

snps_flybase <- str_c(chromosome, 
                      str_c(position, position, sep = "-"), sep = ":")

# save to a .csv file to copy into flybase feature mapper 
write.table(snps_flybase, 
            file = "snps.csv", 
            sep = ",", 
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)

# maybe consider making this whole thing into a function

# copy and paste into flybase feature mapper website ---------------------------
# create a new dataframe from the GFF file

gff <- read.delim("snps_flybase_features.txt", header = FALSE)


# maybe just try it from the GeneAnnotation from the gwas.top.annot files ------

# get the FBgn# for all the snps
min_gwas_top$FBgn <- str_extract(min_gwas_top$GeneAnnotation, "FBgn[[:digit:]]+")

library("AnnotationDbi")
library("org.Dm.eg.db")

# FBgn to gene_symbol
min_gwas_top$gene_symbol <- mapIds(org.Dm.eg.db, 
                                   keys = min_gwas_top$FBgn, 
                                   column = "SYMBOL", 
                                   keytype = "FLYBASE",
                                   multiVals = "first")

# get the 
str_replace(min_gwas_top$GeneAnnotation[1], "^TranscriptAnnot\\[", "")


