#using the Ecological Genomics script to run DESeq

#setwd("~/github/2020_Ecological_Genomics")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)

#write.csv(conds, "\\Users\\lockwoodlab\\Desktop\\RNA Seq\\conds.csv", row.names = TRUE)

countsTable <- read.table("Dm_countsMatrix.txt", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
countsTableRound <- round(countsTable) # Need to round because DESeq wants only integers
head(countsTableRound)

colSums(countsTableRound)
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound), las=3, cex.names=0.5,names.arg = substring(colnames(countsTableRound),1,13)) 
abline(h=mean(colSums(countsTableRound)), col="blue", lwd =2)

z <-lapply(countsTableRound, function(countsTableRound){ length(which(countsTableRound==0))})
str(z)

# 

# Need the colClasses otherwise imports day as numeric which DESeq doesn't like
conds <- read.delim("DM_samples2.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1, colClasses=c('factor', 'factor', 'factor'))
head(conds)

int <- data.frame(intersect(colnames(countsTable), rownames(conds)))



######################################

#dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, design = ~  temp)
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, design = ~ pop + temp + pop:temp)
# "treatment effect" model controling for differences in day and population

dds$temp <- factor(dds$temp, levels = c("25","32", "34", "36"))
dds$pop <- factor(dds$pop, levels = c("BO", "CP", "GH", "GM", "SK", "VT2", "VT8", "VT9", "VT10", "VT12"))

dim(dds)
# 66069    76

#dds <- dds[rowSums(counts(dds)) > 100,]

dim(dds)
# 38340    76  At > 1, More reasonable number, lost about 42%
# 9771   76    At >100


#dds$group <- factor(paste0(dds$pop, dds$temp))
#design(dds) <- ~ group
#dds <- DESeq(dds)
#resultsNames(dds)
#results(dds, contrast=c("pop", "temp"))


dds <- DESeq(dds, modelMatrixType = "standard")

resultsNames(dds)
# [1] "Intercept"        "day_10_vs_0"      "day_5_vs_0"       "climate_HD_vs_CW"
# [5] "treatment_D_vs_C" "treatment_H_vs_C"

res <- results(dds)
str(res)

res <- res[order(res$padj),]
head(res, l=5,5)
# log2 fold change (MLE): treatment H vs C 
# Wald test p-value: treatment H vs C 
# DataFrame with 5 rows and 6 columns
# baseMean    log2FoldChange            lfcSE              stat
# <numeric>         <numeric>        <numeric>         <numeric>
#   MA_47071g0010     1.17999411601304  18.5693308669078 2.02440183914411  9.17274945509766
# MA_35043g0020    0.215040900686398 -29.5556958585507 3.58083855075752 -8.25384765037736
# MA_33964g0010    0.927024448716801  -28.738155543989 3.58091582986403 -8.02536471377498
# MA_10433150g0020 0.839433287586582 -28.4792813092653  3.5796693586627 -7.95584129588566
# MA_10813g0010    0.461606021497004  28.0506299312649 3.58158027900867  7.83191433559852
# pvalue                 padj
# <numeric>            <numeric>
#   MA_47071g0010     4.6110524296846e-20 1.69972614663034e-15
# MA_35043g0020    1.53374848780567e-16 2.82685183787463e-12
# MA_33964g0010    1.01224754712075e-15 1.24378230273216e-11
# MA_10433150g0020 1.77918345440023e-15 1.63960651240253e-11
# MA_10813g0010    4.80496893372401e-15 2.53029664049906e-11

summary(res)

# out of 38340 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 68, 0.18%
# LFC < 0 (down)     : 42, 0.11%
# outliers [1]       : 0, 0%
# low counts [2]     : 1478, 3.9%
# (mean count < 0)

res_treatCD <- results(dds, name="treatment_D_vs_C", alpha=0.05)

res_treatCD <- res_treatCD[order(res_treatCD$padj),]
summary(res_treatCD)
# out of 38340 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 154, 0.4%
# LFC < 0 (down)     : 75, 0.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 4460, 12%
# (mean count < 0)

res_treatCH <- results(dds, name="treatment_H_vs_C", alpha=0.05)

res_treatCH <- res_treatCH[order(res_treatCH$padj),]
summary(res_treatCH)  # Thought it would be the same number as first res above...
# out of 38340 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 61, 0.16%
# LFC < 0 (down)     : 37, 0.097%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)

################  Data visualization
plotMA(res_treatCH, main="DESeq2", ylim=c(-2,2))
abline(h=c(-1,1), col="blue", lwd =2)

# hot & dry effect

res_treatCD <- results(dds, name="treatment_D_vs_C", alpha=0.05)
plotMA(res_treatCD, main="DESeq2", ylim=c(-4,4))  # What shape are we looking for in MA plot?
abline(h=c(-1,1), col="blue", lwd =2)

#####  PCA
###########################
#Thomas making PCA function
x <- counts(dds, normalized = TRUE)

x <- x[apply(x,1,var) != 0, ]
x <- t(x)

x %>%
  mutate(var = )
  filter()

pca <- prcomp(t(x), scale = TRUE)

plot(pca$x[,1], pca$x[,2])
##########################
rld <- rlog(dds, blind=FALSE)
vsd <- vst(dds, blind=FALSE)


data <- plotPCA(vsd, ntop= 29868, intgroup=c("pop","temp"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

data$temp <- factor(data$temp, levels=c("25","32","34", "36"), labels = c("25","32","34", "36"))
data$pop <- factor(data$pop, levels=c("BO","CP","GH", "GM", "SK", "VT2", "VT8", "VT9", "VT10", "VT12"), labels = c("BO","CP","GH", "GM", "SK", "VT2", "VT8", "VT9", "VT10", "VT12"))

data <- data %>%
  mutate(region = case_when(pop %in% c("BO", "CP", "GH", "GM", "SK") ~ "trop",
                            !(pop %in% c("BO", "CP", "GH", "GM", "SK")) ~ "temperate"))


ggplot(data, aes(PC1, PC2, color=pop, shape=temp)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

ggplot(data, aes(PC1, PC2, color=region, shape=temp)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()


#using pcaplot function

PCAPlot(vsd, group=data$region)




#ggplot(data, aes(PC1, PC2, color=climate, shape=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  theme_minimal() + theme(text = element_text(size=15), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=2))
dev.off()

############# the next several lines will make a cluster heatmap of all the samples vs all samples

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$climate, vsd$treatment, vsd$day, sep="-")
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette(rev(brewer.pal(9,"Oranges")))(255) # 255 delimits the number of bins in the color scale

colors <- colorRampPalette(c("purple","yellow"))(255)

library("pheatmap")

pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists, col=colors)

###################### Let's look at individual genes!!
d <-plotCounts(dds, gene="MA_47071g0010", intgroup = (c("treatment","day","climate")), returnData=TRUE)
d

p <-ggplot(d, aes(x=climate, y=count, shape=day, colour = treatment)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("CW","HD"))
p
