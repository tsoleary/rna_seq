# DESeq workflow in R ----------------------------------------------------------

# Load DESeq2 library
library("DESeq2")

# Set the working directory
directory <- "/Users/tsoleary/R/rna_seq/DM6_counts"
setwd(directory)

# Set the prefix for all output file names
prefix <- "Dm_DESeq2"

sampleTable <- read.delim("/Users/tsoleary/R/rna_seq/DM6_counts/target.txt")

colnames(sampleTable) <- c("sampleFiles", "sampleNames", "sampleCondition")

treatments <- as.character(unique(sampleTable$sampleCondition))

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ sampleCondition)

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$sampleCondition,
                                      levels = treatments)

# DESeq2 normalization ---------------------------------------------------------
dds <- DESeq(ddsHTSeq)

# uncomment the particular comparison
res <- res_hot_con <- results(dds, contrast = c("sampleCondition", "HOT","CON"))
# res <- res_cold_con <- results(dds, contrast = c("sampleCondition", "COLD","CON"))
# res <- res_hot_cold <- results(dds, contrast = c("sampleCondition", "HOT","COLD"))
  
# Set the prefix for each output file name
compCond <- "hot_ctrl"
outputPrefix <- paste(prefix, compCond, sep = "_")

# order results by padj value (most significant to least)
# res <- subset(res, padj < 0.05)
res <- res[order(res$padj), ]  # should see df of baseMean, log2Foldchange, stat, pval, padj

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)

names(resdata)[1] <- "gene"

# write results to file for each set of camparisons
setwd("/Users/tsoleary/R/rna_seq/results")
write.csv(resdata, file = paste0(outputPrefix, "_results_with_normalized.csv"), 
          row.names = FALSE)

# only one normalized counts file created --------------------------------------
# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), 
            file = paste0(prefix, "_normalized_counts.txt"), sep = '\t')


