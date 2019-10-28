# gene set enrichment analysis -------------------------------------------------
# https://pt.coursera.org/lecture/statistical-genomics/gene-set-analysis-in-r-7-43-REPHH


library(devtools)
library(Biobase)
library(DESeq2)
library(goseq)

# Set the working directory
directory <- "/Users/tsoleary/R/rna_seq/DM6_counts"
setwd(directory)

# Set the prefix for all output file names
prefix <- "Dm_DESeq2"

# need results(dds) object
# source("/Users/tsoleary/R/rna_seq/scripts/DESeq_r.R")

res <- res_hot_con
# res <- res_cold_con
# res <- res_hot_cold

genes <- as.integer(res$padj < 0.05)
not_na <- !is.na(genes)
names(genes) <- rownames(res)
genes <- genes[not_na]

pwf <- nullp(genes, "dm3", id = "geneSymbol")

GO.genes <- goseq(pwf, "dm3", id = "geneSymbol") 

head(GO.genes)





