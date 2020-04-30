# ------------------------------------------------------------------------------
# Frontiers reviewer comments -- Gene based GWAS set up for Seth
# April 30, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)

# Function to pull gene based annotation ----

# ------------------------------------------------------------------------------
# Function: annotate_gene_dgrp
# Description: Annotate DGRP SNPs to genes
# Inputs: data.frame with a GeneAnnotation column
# Outputs: data.frame with added fbgn, gene, feature, and pos columns

require(tidyverse)

annotate_gene_dgrp <- function(dat) {
  # extract only the bit at the beginning
  x <- stringr::str_extract(dat$GeneAnnotation, 
                            "SiteClass\\[[\\w\\W|]+\\],")
  
  # remove the SiteClass and brackets etc.
  dat$annot <- stringr::str_remove(stringr::str_remove(x, "SiteClass\\["), 
                                   "\\],")
  
  # split up different genes and save fbgn, gene, feature, and pos in new columns
  df <- dat %>%
    tidyr::separate(col = annot, into = c("gene1", "gene2"), sep = ";") %>%
    tidyr::pivot_longer(gene1:gene2, names_to = "lab", values_to = "gene",
                        values_drop_na = TRUE) %>%
    tidyr::separate(col = gene, 
                    into = c("fbgn", "gene", "feature", "pos"), sep = "\\|") 
  
  return(df)
} 
# End function -----------------------------------------------------------------


# Load necessary files ---------------------------------------------------------

# Load gwas.all.assoc files
ctmin <- read_delim(here::here("DRGP_GWAS/CTmin/gwas.all.assoc"),
                    delim = " ") %>%
  filter(AvgMixedPval <= 0.01)
ctmax <- read_delim(here::here("DRGP_GWAS/CTmax/gwas.all.assoc"),
                    delim = " ") %>%
  filter(AvgMixedPval <= 0.01)


# Load the dgrp.fb557.annot.txt file
snp_annot <- read_delim(here::here("DRGP_GWAS/dgrp.fb557.annot.txt"),
                        delim = "\t", 
                        col_names = c("ID", 
                                      "Allele", 
                                      "GeneAnnotation", 
                                      "RegulatoryAnnotation"))

# Add annotation to the gwas files
ctmin <- ctmin %>%
  left_join(snp_annot)
ctmax <- ctmax %>%
  left_join(snp_annot)

# Save these data frames as a tsv
write_tsv(ctmin, here::here("DRGP_GWAS/ctmin_01_annot.tsv"))
write_tsv(ctmax, here::here("DRGP_GWAS/ctmax_01_annot.tsv"))

# Parse the gene annotations
ctmin <- annotate_gene_dgrp(ctmin) %>%
  select(ID, MinorAllele, MajorAllele, MAF, 
         AvgMixedPval, fbgn, gene, feature, pos)

ctmax <- annotate_gene_dgrp(ctmax) %>%
  select(ID, MinorAllele, MajorAllele, MAF, 
         AvgMixedPval, fbgn, gene, feature, pos)


# Save these data frames with parsed gene annotation as a tsv
write_tsv(ctmin, here::here("DRGP_GWAS/ctmin_01_annot_parsed.tsv"))
write_tsv(ctmax, here::here("DRGP_GWAS/ctmax_01_annot_parsed.tsv"))

