# DESeq workflow in R ----------------------------------------------------------

# Load DESeq2 library
library("DESeq2")

# Set the working directory
directory <- "/Users/tsoleary/R/rna_seq/DM6_counts"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "Dm_DESeq2"

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
res_all <- results(dds)

################################################################################
################################################################################
res_hot_con <- results(dds, contrast=c("condition","HOT","CON"))
res_cold_con <- results(dds, contrast=c("condition","COLD","CON"))
# res_hot_cold <- results(dds, contrast=c("condition","HOT","COLD"))
################################################################################
################################################################################

# order results by padj value (most significant to least)
res <- subset(res, padj < 0.05)
res <- res[order(res$padj), ] 
     # should see DataFrame of baseMean, log2Foldchange, stat, pval, padj

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)

names(resdata)[1] <- "gene"

write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), 
            file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
# mcols(res, use.names = T)
# write.csv(as.data.frame(mcols(res, use.name = T)), 
#           file = paste0(outputPrefix, "-test-conditions.csv"))


# Outlier handling -------------------------------------------------------------
# replacing outlier value with estimated value as predicted by distrubution 
# using "trimmed mean" approach. recommended if you have several replicates per 
# treatment. DESeq2 will automatically do this if you have 7 or more replicates

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
# tab <- table(initial = results(dds)$padj < 0.05,
#              cleaned = results(ddsClean)$padj < 0.05)
# addmargins(tab)
# write.csv(as.data.frame(tab), 
#           file = paste0(outputPrefix, "-replaceoutliers.csv"))
# 
# resClean <- results(ddsClean)
# resClean <- subset(res, padj < 0.05)
# resClean <- resClean[order(resClean$padj),]
# write.csv(as.data.frame(resClean), 
#           file = paste0(outputPrefix, "-replaceoutliers-results.csv"))


# MA plot of RNAseq data for entire dataset ------------------------------------
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
# plotMA(dds, ylim=c(-8,8), main = "RNAseq experiment")
# dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
# dev.off()


# MA plot of RNAseq data for entire dataset ------------------------------------
# https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html

library(ggpubr)

# Add rectangle around labels
ggmaplot(res_all, main = expression("Group 1" %->% "Group 2"),
         fdr = 0.05, fc = 2, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res_all$gene),
         legend = "top", top = 20,
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())

dev.copy(png, paste0(outputPrefix, "-MAplot_ggpubr_analysis.png"))
dev.off()


# transform raw counts into normalized values ----------------------------------
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.

rld <- rlogTransformation(dds, blind = TRUE)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# save normalized values
# write.table(as.data.frame(assay(rld),
#                           file = paste0(outputPrefix, 
#                           "-rlog-transformed-counts.txt"), 
#                           sep = '\t'))
# 
# write.table(as.data.frame(assay(vsd), file = paste0(outputPrefix, 
#                                       "-vst-transformed-counts.txt"), 
#                                        sep = '\t'))

# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
# par(mai = ifelse(1:4 <= 2, par('mai'),0))
# px <- counts(dds)[,1] / sizeFactors(dds)[1]
# ord <- order(px)
# ord <- ord[px[ord] < 150]
# ord <- ord[seq(1,length(ord),length=50)]
# last <- ord[length(ord)]
# vstcol <- c('blue','black')
# matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, 
#         col=vstcol, xlab = 'n', ylab = 'f(n)')
# legend('bottomright',legend=c(expression('variance stabilizing transformation'), 
#                               expression(log[2](n/s[1]))), fill=vstcol)
# dev.copy(png,paste0(outputPrefix, "-variance_stabilizing.png"))
# dev.off()

# clustering analysis ----------------------------------------------------------
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
library("RColorBrewer")
library("gplots")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames, sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(10, 10))
dev.copy(png, paste0(outputPrefix, "-clustering.png"))
dev.off()

# Principal components plot analysis -------------------------------------------
# shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

# set condition
condition <- treatments
scores <- data.frame(pc$x, condition)

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.9,.2),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))

ggsave(pcaplot,file=paste0(outputPrefix, "-pca-plot-ggplot2.pdf"))


# scatter plot of rlog transformations between Sample conditions ---------------
# nice way to compare control and experimental samples
head(assay(rld))
plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3,
     main='Log2 transformed')

par(mfrow=c(2,2))
plot(assay(rld)[,c(2,1)],col='#00000020',pch=20,cex=0.5)
plot(assay(rld)[,c(3,1)],col='#00000020',pch=20,cex=0.5)
plot(assay(rld)[,c(4,1)],col='#00000020',pch=20,cex=0.5)
plot(assay(rld)[,c(5,1)],col='#00000020',pch=20,cex=0.5)

par(mfrow=c(3,3))





# heatmap of data
library("RColorBrewer")
library("gplots")
# 1000 top expressed genes with heatmap.2
select <- order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:500]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="500 Top Expressed Genes Heatmap")


dev.copy(png, paste0(outputPrefix, "-HEATMAP.png"))
dev.off()

library(pheatmap)

pheatmap(assay(vsd)[select,], color = my_palette,
         show_rownames = FALSE, scale = "row",
          main="1000 Top Expressed Genes Heatmap")

dev.copy(png, paste0(outputPrefix, "-PHEATMAP.png"))
dev.off()

