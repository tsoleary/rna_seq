# Principle Component Analysis -------------------------------------------------

# read in the normalized counts .csv file (obtained from PIVOT/DESeq norm)
norm_counts <- read.csv("/Users/tsoleary/R/rna_seq/data/dm_counts_norm.csv")

# matrix of values
data.matrix <- as.matrix(norm_counts[, -1])

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

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample = rownames(pca$x),
                       X = pca$x[,1],
                       Y = pca$x[,2])

pca.data$group <- gsub("_[[:alnum:]]", "",  pca.data$Sample)

# plot with sample labels
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA")

# plot colored by group
ggplot(data = pca.data, aes(x = X, y = Y, fill = group)) +
  geom_point(pch = 21, size =5) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA")


# specific genes that are contributing the most to pc1
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
names(gene_scores) <- norm_counts$X
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes, 1] ## show the scores (and +/- sign)

# lets see what everything looks like without con_3 and con_4 ------------------

norm_counts <- select(norm_counts, -c("con_3", "con_4"))

# re-run above PCA code
