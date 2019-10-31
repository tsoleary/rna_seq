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

