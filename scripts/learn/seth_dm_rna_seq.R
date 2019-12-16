# raw read count files from Seth

library(dplyr)

samples <- read.delim("/Users/tsoleary/R/rna_seq/target.txt")

setwd("/Users/tsoleary/R/rna_seq/DM6_counts")

list.files()

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


write.csv(x, "dm_rna_seq_counts.csv", row.names = FALSE)

library(PIVOT)
pivot()

# Try DEseq2 normalization by hand later

# Principle Component Analysis -------------------------------------------------
setwd("/Users/tsoleary/R/rna_seq")
list.files()
norm_counts <- read.csv("dm_counts_norm.csv")

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


# differentially expressed genes -----------------------------------------------

genes <- read.csv("dm_sig_genes.csv")

colnames(genes)[1] <- "gene"

plot(genes$padj, genes$log2FoldChange)

sort(abs(genes$log2FoldChange), decreasing = TRUE)[1:59]

fun_genes <- genes[which(abs(genes$log2FoldChange) >= 13.31385), ]

list_genes <- c()

for (i in 1:length(fun_genes$gene)) {
  
  x <- which(norm_counts$X == as.character(fun_genes$gene)[i])

  list_genes[i] <- x
  
}

asdf <- norm_counts[list_genes, ]

jkl <- as.matrix(asdf[asdf$X == "Sgs7", 2:10])

barplot(jkl)


# chromosomal locations ------------------
library(org.Dm.eg.db)



x <- org.Dm.egSYMBOL
# Get the entrez gene identifiers that are mapped to chromosome locations
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
# Get the CHRLOC for the first five genes
xx[1:5]
# Get the first one
xx[[1]]

columns(org.Dm.eg.db)

ann <- biomaRt::select(org.Dm.eg.db, keys = as.character(norm_counts$X), 
                       columns = c("ENTREZID", "SYMBOL", "GENENAME"), keytype = )

norm_counts$X
