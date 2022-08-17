# ------------------------------------------------------------------------------
# Brent Annotate Por's data
# March 18, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Source functions
source("gene_annot_funcs.R")

# Load data --------------------------------------------------------------------
# This data is just a 1000 rows of a table taken from a Por's earlier
# results, as an example of what format the functions expect of the chromosomal 
# positions. 
dat <- readRDS("reprex_snp_dat.rds")

# Load most current gtf genome annotation with corresponding release -----------
# ftp://ftp.flybase.org/genomes/dmel/current/gtf/
gtf_df <- read_clean_gtf("dmel-all-r6.39.gtf", bp_flanking = 1000)

# Annotate the genes and features associated with each snp ---------------------
# This function adds a gene_assoc colummn with gene symbol and 
# gene features (e.g. exon) separated by a semicolon < ; > and any additional 
# associated genes separated by a pipe symbol < | >
dat <- gene_assoc_snp(dat, gtf_df)



