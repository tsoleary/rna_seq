# Cahan Lab DEG ----------------------------------------------------------------

require(dplyr)
require(DESeq)
require(ggpubr)
require(ggplot2)
require(VennDiagram)

directory_counts <- here::here("cahan_counts")
directory_results <- here::here("cahan_results")
directory_plots <- here::here("cahan_plots")

prefix <- "Dm_cahan_deg"

# Import in .counts files ------------------------------------------------------

# working directory containg the .counts files
setwd(directory_counts)

# sample, group, filename information dataframe
metadata <- read.csv("samples.txt")

# read in each .cnt and compile a dataframe of read counts for each sample
count_df <- data.frame()

for (i in 1:length(as.character(metadata$file))){
  t <- read.table(as.character(metadata$file)[i], header = TRUE)
  t <- t[, c(1, length(t))]
  colnames(t) <- c("gene", as.character(metadata$label)[i])
  if (i == 1){
    count_df <- t 
  } else{
    count_df <- full_join(count_df, t, by = "gene")
  }
}

# create .csv file to upload to PIVOT for DESeq normalization
write.csv(count_df, "dm_seq_counts.csv", row.names = FALSE)

counts_table <- as.matrix(count_df[,2:length(count_df)])
rownames(counts_table) <- count_df$gene

# DESeq workflow in R ----------------------------------------------------------

# Set the prefix for each output file name
prefix <- "Dm_DESeq2"

colnames(metadata) <- c("label", "file", "condition")

treatments <- as.character(unique(metadata$condition))

ddsHTSeq <- DESeqDataSetFromMatrix(countData = counts_table,
                                   colData = metadata,
                                   design = ~ condition)

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

# DESeq2 normalization ---------------------------------------------------------

dds <- DESeq(ddsHTSeq)

# set the ctrl group as the control to be compared against
dds$condition <- relevel(dds$condition, ref = "ctrl")

# get the results of the differential expression
res_all <- results(dds)
res_hot <- results(dds, contrast = c("condition", "hot", "ctrl"))
res_cold <- results(dds, contrast = c("condition", "cold", "ctrl"))


# set data results and normalized reads to df
resdata_hot <- merge(as.data.frame(res_hot), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)
names(resdata_hot)[1] <- "gene"

resdata_cold <- merge(as.data.frame(res_cold), 
                      as.data.frame(counts(dds, normalized = TRUE)), 
                      by = "row.names", 
                      sort = FALSE)
names(resdata_cold)[1] <- "gene"

# write results to file for each set of camparisons
setwd(directory_results)

write.csv(resdata_hot, 
          file = paste0(prefix, "_hot", "_results_with_norm.csv"), 
          row.names = FALSE)
write.csv(resdata_cold, 
          file = paste0(prefix, "_cold", "_results_with_norm.csv"), 
          row.names = FALSE)

# write normalized results 
write.table(as.data.frame(counts(dds, normalized = TRUE)), 
            file = paste0(prefix, "_norm_counts.txt"), sep = '\t')

# MA plot of RNAseq data for entire dataset ------------------------------------

# groups being compared as how it was saved in the result folder
comp_cond <- "hot"
output_prefix <- paste(prefix, comp_cond, sep = "_")

# load result file
setwd(directory_results)
csv_result_file <- list.files()[grepl(comp_cond, list.files())]
res <- read.csv(csv_result_file)

# create ma plot
ggmaplot(res, main = comp_cond,
         fdr = 0.05, fc = 2, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(res$gene),
         legend = "top", top = 20,
         select.top.method = "padj",
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) 

setwd(directory_plots)
dev.copy(png, paste0(prefix, "_MAplot_ggpubr_analysis.png"))
dev.off()

# Principle Componenet Analysis ------------------------------------------------

# load norm_counts file
setwd(directory_results)
csv_result_file <- list.files()[grepl("norm_counts", list.files())]
norm_counts <- read.table(csv_result_file)

# matrix of values
data.matrix <- as.matrix(norm_counts)
data.matrix <- data.matrix[apply(data.matrix, 1, var) != 0, ]
pca <- prcomp(t(data.matrix), scale = TRUE) 

## plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main = "Scree Plot", 
        xlab = "Principal Component", 
        ylab = "Percent Variation")

# now make a fancy looking plot that shows the PCs and the variation:
pca.data <- data.frame(Sample = rownames(pca$x),
                       X = pca$x[,1],
                       Y = pca$x[,2])

pca.data$group <- gsub("_[[:alnum:]]", "",  pca.data$Sample)

# plot with sample labels
ggplot(data = pca.data, aes(x = X, y = Y, label = Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA")

setwd(directory_plots)
dev.copy(png, paste0(prefix, "_pca_text.png"))
dev.off()

# plot colored by group
ggplot(data = pca.data, aes(x = X, y = Y, fill = group)) +
  geom_point(pch = 21, size = 5) +
  scale_fill_manual(values=c("skyblue", "#999999", "firebrick1")) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA") 

setwd(directory_plots)
dev.copy(png, paste0(prefix, "_pca_colored.png"))
dev.off()



# Venn diagram of DESeq results ------------------------------------------------

hot_deg <- row.names(subset(res_hot, res_hot$padj <= 0.05))

cold_deg <- row.names(subset(res_cold, res_cold$padj <= 0.05))

setwd(directory_plots)

venn.diagram(x = list(HOT = hot_deg, COLD = cold_deg),
             filename = paste0(prefix, "_venn_diagram.tiff"),
             lwd = 4, 
             fill = c("red", "blue"), 
             alpha = 0.5,
             label.col = "white", 
             cex = 3,
             fontfamily = "serif", 
             fontface = "bold",
             cat.col = c("red", "blue"), 
             cat.cex = 3, 
             cat.fontfamily = "serif",
             cat.fontface = "bold", 
             cat.dist = c(0.03, 0.03), 
             cat.pos = c(200, 180))

# subset to only genes that are significantly upregulated L2FC > 1.5
hot_deg_up <- row.names(subset(res_hot, res_hot$padj <= 0.05 & 
                               res_hot$log2FoldChange > 1.5))

cold_deg_up <- row.names(subset(res_cold, res_cold$padj <= 0.05 &
                                res_cold$log2FoldChange > 1.5))

setwd(directory_plots)

venn.diagram(x = list(HOT = hot_deg_up, COLD = cold_deg_up),
             filename = paste0(prefix, "_up_reg_venn_diagram.tiff"),
             lwd = 4, 
             fill = c("red", "blue"), 
             alpha = 0.5,
             label.col = "white", 
             cex = 3,
             fontfamily = "serif", 
             fontface = "bold",
             cat.col = c("red", "blue"), 
             cat.cex = 3, 
             cat.fontfamily = "serif",
             cat.fontface = "bold", 
             cat.dist = c(0.03, 0.03), 
             cat.pos = c(200, 180))

