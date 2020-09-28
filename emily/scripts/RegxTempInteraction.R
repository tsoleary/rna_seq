library(tidyverse)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(lmtest)
#library(DEGreport)

#For non-arranged samples order
#counts <- round(read.table("Dm_countsMatrix.txt", header = TRUE))

conds <- read.delim("DM_samples2.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1, colClasses=c('factor', 'factor', 'factor'))


# Rearrange the order for the heat map based on region and temperature
# conds <- conds %>%
#   arrange(temp, region) %>%
#   mutate(sample = paste(group = paste(pop, temp, rep, sep = "_")))

# conds <- conds %>%
#   arrange(region, temp) %>%
#   mutate(sample = paste(group = paste(pop, temp, rep, sep = "_")))


# 
# # Counts matrix column order
# counts <- round(read.table("Dm_countsMatrixTO.txt", header = TRUE))
# counts <- counts %>%
#   select(conds$sample)



# Filter based on transcript presence in at least per_samp% of samples
# per_samp <- 0.6

# # creates a vector of transcripts to keep
# keep_genes <- counts %>%
#   rownames_to_column("gene") %>%
#   pivot_longer(contains("_"),
#                names_to = "sample",
#                values_to = "counts") %>%
#   group_by(gene) %>%
#   filter(counts != 0) %>%
#   count() %>%
#   filter(n > per_samp * 110)

# filters based on keep_genes
# counts <- counts %>%
#   rownames_to_column("gene") %>%
#   filter(gene %in% keep_genes$gene) %>%
#   column_to_rownames("gene")


#dds <- DESeqDataSetFromMatrix(countData = counts, colData = conds, design = ~  temp + region + temp:region)
#dim(dds)
#dds <- DESeq(dds)

#using LTR

# design(dds) <- ~ temp + region + temp:region
# design(dds) <- ~  region 
# design(dds) <- ~ temp
# design(dds) <- ~ temp + region

#ddslrt <- DESeq(dds, test = "LRT", reduced = ~ temp + region + temp:region)

#dds <- nbinomLRT(dds, reduced = ~ 1)

#resultsNames(dds)
# [1] "Intercept"                    "temp_32_vs_25"               
# [3] "temp_34_vs_25"                "temp_36_vs_25"               
# [5] "region_tropical_vs_temperate" "temp32.regiontropical"       
# [7] "temp34.regiontropical"        "temp36.regiontropical" 

#LRT temp + region +temp:region
dds <- DESeqDataSetFromMatrix(countData = counts, colData = conds, design = ~  temp + region + temp:region)
ddslrtint <- DESeq(dds, test = "LRT",  reduced = ~temp + region)
resultsNames(ddslrtint)
temp32.regiontropical <- results(ddslrtint, name="temp32.regiontropical", alpha=0.05)
temp34.regiontropical <- results(ddslrtint, name="temp34.regiontropical", alpha=0.05)
temp36.regiontropical <- results(ddslrtint, name="temp36.regiontropical", alpha=0.05)
summary(temp32.regiontropical)
summary(temp34.regiontropical)
summary(temp36.regiontropical)

#lrt temp + region
dds <- DESeqDataSetFromMatrix(countData = counts, colData = conds, design = ~  temp + region)
#ddslrttempreg <- DESeq(dds, test = "LRT",  reduced = ~temp)
ddslrttempregT <- DESeq(dds, test = "LRT",  reduced = ~region)
#ddslrttempreg <- DESeq(dds, test = "LRT",  reduced = ~1)
resultsNames(ddslrttempreg)
resultsNames(ddslrttempregT)
temp_32_vs_25 <- results(ddslrttempregT, name="temp_32_vs_25", alpha=0.05)
#write.csv(temp_32_vs_25, "/Users/lockwoodlab/Desktop/RNA Seq//temp_32_vs_25_tempreg.csv")
summary(temp_32_vs_25)

temp_34_vs_25 <- results(ddslrttempreg, name="temp_34_vs_25", alpha=0.05)
rntemp_34_vs_25 <- rownames(temp_34_vs_25)
temp_36_vs_25 <- results(ddslrttempreg, name="temp_36_vs_25", alpha=0.05)
#get p<0.05 list of transcripts 
temp_36_vs_25_genes <- unique(rownames(temp_36_vs_25)[temp_36_vs_25$padj<0.05])[!is.na(unique(rownames(temp_36_vs_25)[temp_36_vs_25$padj<0.05]))]
#write.csv(temp_36_vs_25_genes, "/Users/lockwoodlab/Desktop/RNA Seq//temp_36_vs_25_tempreggenes.csv")

region_tropical_vs_temperate <- results(ddslrttempreg, name="region_tropical_vs_temperate", alpha=0.05)
region_tropical_vs_temperate_genes <- unique(rownames(region_tropical_vs_temperate)[region_tropical_vs_temperate$padj<0.05])[!is.na(unique(rownames(region_tropical_vs_temperate)[region_tropical_vs_temperate$padj<0.05]))]

summary(temp_32_vs_25)
summary(temp_34_vs_25)

summary(region_tropical_vs_temperate)

results_34_32 <- results(ddslrttempreg, contrast=c("temp", "34", "32"), alpha=0.05)
summary(results_34_32)
head(results_34_32[order(results_34_32$padj), ])

results_36_32 <- results(ddslrttempreg, contrast=c("temp", "36", "32"), alpha=0.05)
summary(results_36_32)
head(results_36_32[order(results_36_32$padj), ])

results_36_34 <- results(ddslrttempreg, contrast=c("temp", "36", "34"), alpha=0.05)
summary(results_36_34)
head(results_36_34[order(results_36_34$padj), ])





#lrt temp
dds <- DESeqDataSetFromMatrix(countData = counts, colData = conds, design = ~  temp )
ddslrttemp <- DESeq(dds, test = "LRT",  reduced = ~1)
resultsNames(ddslrttemp)

temp_32_vs_25 <- results(ddslrttemp, name="temp_32_vs_25", alpha=0.05)
temp_34_vs_25 <- results(ddslrttemp, name="temp_34_vs_25", alpha=0.05)
temp_36_vs_25 <- results(ddslrttemp, name="temp_36_vs_25", alpha=0.05)
summary(temp_32_vs_25)
summary(temp_34_vs_25)
summary(temp_36_vs_25)

write.csv(temp_32_vs_25, "/Users/lockwoodlab/Desktop/RNA Seq//temp_32_vs_25_TEMP.csv")

head(temp_32_vs_25[order(temp_32_vs_25$padj), ])
head(temp_34_vs_25[order(temp_34_vs_25$padj), ])



temp_32_vs_25_df <- rownames_to_column(as.data.frame(temp_32_vs_25), 
                                       var = "transcript") %>%
  filter(padj < 0.05)
temp_34_vs_25_df <- rownames_to_column(as.data.frame(temp_34_vs_25), 
                                       var = "transcript") %>%
  filter(padj < 0.05)
temp_36_vs_25_df <- rownames_to_column(as.data.frame(temp_36_vs_25), 
                                       var = "transcript") %>%
  filter(padj < 0.05)

#bind_rows(temp_32_vs_25_df, temp_34_vs_25_df, temp_36_vs_25_df, .id = "id")

join_temp <- full_join(temp_32_vs_25_df, temp_34_vs_25_df, 
                       by = "transcript", 
                       suffix = c(".32_25", ".34_25"))

joined_temp_results <- full_join(join_temp, temp_36_vs_25_df, 
                                 by = "transcript", 
                                 suffix = c("", ".36_25"))


xjoin <- joined_temp_results %>%
  mutate(shared = case_when(
    log2FoldChange.32_25 > 0 & log2FoldChange.34_25 > 0 & log2FoldChange > 0 ~ "up_all",
    log2FoldChange.32_25 < 0 & log2FoldChange.34_25 < 0 & log2FoldChange < 0 ~ "down_all",
    log2FoldChange.32_25 > 0 & log2FoldChange.34_25 < 0 & log2FoldChange < 0 ~ "up_32_down_rest",
    log2FoldChange.32_25 < 0 & log2FoldChange.34_25 > 0 & log2FoldChange < 0 ~ "up_34_down_rest",
    log2FoldChange.32_25 > 0 & log2FoldChange.34_25 > 0 & log2FoldChange < 0 ~ "up_36_down_rest",
    log2FoldChange.32_25 < 0 & log2FoldChange.34_25 > 0 & log2FoldChange > 0 ~ "down_32_up_rest",
    log2FoldChange.32_25 > 0 & log2FoldChange.34_25 < 0 & log2FoldChange > 0 ~ "down_34_up_rest",
    log2FoldChange.32_25 < 0 & log2FoldChange.34_25 < 0 & log2FoldChange > 0 ~ "down_32_34_up_36"))
    
    
    
    # log2FoldChange.32_25 > 0 & log2FoldChange.34_25 > 0 & log2FoldChange < 0 ~ "up_32_34_down_36",
    # 
    # log2FoldChange.32_25 > 0 & log2FoldChange.34_25 < 0 & log2FoldChange > 0 ~ "up_32_36_down_34",
    # log2FoldChange.32_25 < 0 & log2FoldChange.34_25 > 0 & log2FoldChange < 0 ~ "down_32_36_up_34"))

xjoin %>%
  group_by(shared) %>%
  tally()

y <- xjoin %>%
  filter(shared == "down_32_34_up_36")

Upall_temp <- y$transcript
downall_temp <- y$transcript
up_32_down_rest_temp <- y$transcript
down32 <- y$transcript
up34 <- y$transcript
up36 <- y$transcript
down36 <- y$transcript


#TEMP_temp_32_vs_25_genes <- unique(rownames(temp_32_vs_25)[temp_32_vs_25$padj<0.05])[!is.na(unique(rownames(temp_32_vs_25)[temp_32_vs_25$padj<0.05]))]


results_34_32 <- results(ddslrttemp, contrast=c("temp", "34", "32"), alpha=0.05)
summary(results_34_32)
head(results_34_32[order(results_34_32$padj), ])

results_36_32 <- results(ddslrttemp, contrast=c("temp", "36", "32"), alpha=0.05)
summary(results_36_32)
head(results_36_32[order(results_36_32$padj), ])

results_36_34 <- results(ddslrttemp, contrast=c("temp", "36", "34"), alpha=0.05)
summary(results_36_34)
head(results_36_34[order(results_36_34$padj), ])




#lrt region
keep <- rowMeans(counts(dds)) >= 5
dds <- dds[keep,]
#dds <- dds[rowSums(counts(dds)) > 220,]
dds <- DESeqDataSetFromMatrix(countData = counts, colData = conds, design = ~  region)
ddslrtreg <- DESeq(dds, test = "LRT",  reduced = ~1)
resultsNames(ddslrtreg)
region_tropical_vs_temperateA <- results(ddslrtreg, name="region_tropical_vs_temperate", alpha=0.05)
summary(region_tropical_vs_temperateA)
write.csv(region_tropical_vs_temperateA, "/Users/lockwoodlab/Desktop/RNA Seq//tropvtempREGION.csv")


tropsig <-read.csv("UpTropical_region.csv", header=TRUE)
tropsigT<- tropsig$transcript

tempsig <-read.csv("UpTemperate_region.csv", header=TRUE)
tempsigT<- tempsig$transcript

# Reg_temp_32_vs_25_genes <- unique(rownames(region_tropical_vs_temperate)[region_tropical_vs_temperate$padj<0.05])[!is.na(unique(rownames(region_tropical_vs_temperate)[region_tropical_vs_temperate$padj<0.05]))]


#"join" lists of DE transcripts to check if overlapping
#Region
x <- rownames_to_column(as.data.frame(region_tropical_vs_temperate), var = "transcript") %>%
  filter(padj < 0.05)
xA <- rownames_to_column(as.data.frame(region_tropical_vs_temperateA), var = "transcript") %>%
  filter(padj < 0.05)

x_join <- full_join(x, xA, by = "transcript")

sum(!is.na(x_join$padj.regwtemp))

#Temperature
xT <- rownames_to_column(as.data.frame(temp_32_vs_25), var = "transcript") %>%
  filter(padj < 0.05)
xTA <- rownames_to_column(as.data.frame(temp_32_vs_25A), var = "transcript") %>%
  filter(padj < 0.05)

x_join <- full_join(xT, xTA, by = "transcript", suffix = c(".tempwreg", "temp"))

length(xT$transcript) 
length(xTA$transcript) - length(intersect(xT$transcript, xTA$transcript))

length(intersect(xT$transcript, xTA$transcript))


sum(!is.na(x_join$padj.temp))



#lrt for temp + pop +temp:pop
dds <- DESeqDataSetFromMatrix(countData = counts, colData = conds, design = ~  temp + pop + temp:pop)
ddslrtintpop <- DESeq(dds, test = "LRT",  reduced = ~temp + pop)
resultsNames(ddslrtintpop)

res_pop_VT10_vs_BO <- results(ddslrtintpop, name="pop_VT10_vs_BO", alpha=0.05)
summary(res_pop_VT10_vs_BO)
head(res_pop_VT10_vs_BO[order(res_pop_VT10_vs_BO$padj), ])
write.csv(res_pop_VT10_vs_BO, "/Users/lockwoodlab/Desktop/RNA Seq//res_pop_VT10_vs_BO.csv")

res_temp34.popVT9 <- results(ddslrtintpop, name="temp34.popVT9", alpha=0.05)
summary(res_temp34.popVT9)
head(res_temp34.popVT9[order(res_temp34.popVT9$padj), ])

res_temp36.popVT9 <- results(ddslrtintpop, name="temp36.popVT9", alpha=0.05)
summary(res_temp36.popVT9)
head(res_temp36.popVT9[order(res_temp36.popVT9$padj), ])


res_temp34.popVT9_genes <- unique(rownames(res_temp34.popVT9)[res_temp34.popVT9$padj<0.05])[!is.na(unique(rownames(res_temp34.popVT9)[res_temp34.popVT9$padj<0.05]))]


#lrt for temp + pop 
dds <- DESeqDataSetFromMatrix(countData = counts, colData = conds, design = ~  temp + pop)
ddslrtemppop <- DESeq(dds, test = "LRT",  reduced = ~temp)
resultsNames(ddslrtemppop)

TPres_pop_VT10_vs_BO <- results(ddslrtemppop, name="pop_VT10_vs_BO", alpha=0.05)
summary(TPres_pop_VT10_vs_BO)
head(TPres_pop_VT10_vs_BO[order(TPres_pop_VT10_vs_BO$padj), ])
write.csv(TPres_pop_VT10_vs_BO, "/Users/lockwoodlab/Desktop/RNA Seq//TPres_pop_VT10_vs_BO.csv")

res_pop_VT9_vs_BO <- results(ddslrtemppop, name="pop_VT9_vs_BO", alpha=0.05)
summary(res_pop_VT9_vs_BO)
head(res_pop_VT9_vs_BO[order(res_pop_VT9_vs_BO$padj), ])

res_pop_CP_vs_BO <- results(ddslrtemppop, name="pop_CP_vs_BO", alpha=0.05)
summary(res_pop_CP_vs_BO)
head(res_pop_CP_vs_BO[order(res_pop_CP_vs_BO$padj), ])

res_pop_CP_vs_BO_genes <- unique(rownames(res_pop_CP_vs_BO)[res_pop_CP_vs_BO$padj<0.05])[!is.na(unique(rownames(res_pop_CP_vs_BO)[res_pop_CP_vs_BO$padj<0.05]))]


#lrt for pop
dds <- DESeqDataSetFromMatrix(countData = counts, colData = conds, design = ~  pop)
ddslrtpop <- DESeq(dds, test = "LRT",  reduced = ~1)
resultsNames(ddslrtpop)

POPres_pop_VT10_vs_BO <- results(ddslrtpop, name="pop_VT10_vs_BO", alpha=0.05)
summary(POPres_pop_VT10_vs_BO)
head(POPres_pop_VT10_vs_BO[order(POPres_pop_VT10_vs_BO$padj), ])
write.csv(POPres_pop_VT10_vs_BO, "/Users/lockwoodlab/Desktop/RNA Seq//POPres_pop_VT10_vs_BO.csv")





###############
temp_32_vs_25 <- results(ddslrt, name="temp_32_vs_25", alpha=0.05)
#temp_32_vs_25_1 <- results(dds, name="temp_32_vs_25", alpha=0.001)
temp_34_vs_25 <- results(ddslrt, name="temp_34_vs_25", alpha=0.05)
temp_36_vs_25 <- results(dds, name="temp_36_vs_25", alpha=0.05)
region_tropical_vs_temperate <- results(ddslrt, name="region_tropical_vs_temperate", alpha=0.05)
temp32.regiontropical <- results(ddslrt, name="temp32.regiontropical", alpha=0.05)
temp34.regiontropical <- results(ddslrt, name="temp34.regiontropical", alpha=0.05)
temp36.regiontropical <- results(ddslrt, name="temp36.regiontropical", alpha=0.05)
#write.csv(temp36.regiontropical, "/Users/lockwoodlab/Desktop/RNA Seq//res_int36.csv")





summary(temp_32_vs_25)
head(temp_32_vs_25[order(temp_32_vs_25$padj), ])
summary(temp32.regiontropical)
head(temp_34_vs_25[order(temp_34_vs_25$padj), ])

summary(region_tropical_vs_temperate)
head(region_tropical_vs_temperate[order(region_tropical_vs_temperate$padj), ])



matrix(resultsNames(dds)) #gives same as above but as matrix :) 
# [,1]                          
# [1,] "Intercept"                   
# [2,] "temp_32_vs_25"               
# [3,] "temp_34_vs_25"               
# [4,] "temp_36_vs_25"               
# [5,] "region_tropical_vs_temperate"
# [6,] "temp32.regiontropical"       
# [7,] "temp34.regiontropical"       
# [8,] "temp36.regiontropical" 

#defining pairwise comparisons from contrasts
results_region <- results(ddslrt, contrast=c("region", "tropical", "temperate"))
summary(results_region)
head(results_region[order(results_region$padj), ])

results_34_36 <- results(dds, contrast=c("temp", "34", "36"))
summary(results_34_36)
head(results_34_36[order(results_34_36$padj), ])

#lrt
results_25_32 <- results(ddslrt, contrast=c("temp", "32", "25"), alpha=0.05)
summary(results_25_32)
head(results_25_32[order(results_25_32$padj), ])

results_int_36 <- results(dds, contrast=c("temp:region", "36.tropical", "36.temperate"))

#MA plots

#32°C vs.  25°C
ggmaplot(temp_32_vs_25, main = "32°C vs.  25°C",
         fdr = 0.05, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(temp_32_vs_25$gene),
         legend = "top", top = 0,
         select.top.method = "padj",
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = 12))

summary(temp_32_vs_25)
head(temp_32_vs_25[order(temp_32_vs_25$padj), ])

#34°C vs.  25°C
ggmaplot(temp_34_vs_25, main = "34°C vs.  25°C",
         fdr = 0.05, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(temp_34_vs_25$gene),
         legend = "top", top = 0,
         select.top.method = "padj",
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = 12))

summary(temp_34_vs_25)
head(temp_34_vs_25[order(temp_34_vs_25$padj), ])

#36°C vs.  25°C
ggmaplot(temp_36_vs_25, main = "36°C vs.  25°C",
         fdr = 0.05, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(temp_36_vs_25$gene),
         legend = "top", top = 0,
         select.top.method = "padj",
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = 12))

summary(temp_36_vs_25)
head(temp_36_vs_25[order(temp_36_vs_25$padj), ])

#Tropical vs. Temperate
ggmaplot(region_tropical_vs_temperate, main = "Tropical vs. Temperate",
         fdr = 0.05, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(region_tropical_vs_temperate$gene),
         legend = "top", top = 0,
         select.top.method = "padj",
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = 12))

summary(region_tropical_vs_temperate)
head(region_tropical_vs_temperate[order(region_tropical_vs_temperate$padj), ])

#Tropical vs. Temperate Interaction @ 32°C
ggmaplot(temp32.regiontropical, main = "Tropical vs. Temperate @ 32°C",
         fdr = 0.05, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(temp32.regiontropical$gene),
         legend = "top", top = 0,
         select.top.method = "padj",
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = 12))

summary(temp32.regiontropical)
head(temp32.regiontropical[order(temp32.regiontropical$padj), ])

#Tropical vs. Temperate Interaction @ 34°C
ggmaplot(temp34.regiontropical, main = "Tropical vs. Temperate @ 34°C",
         fdr = 0.05, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(temp34.regiontropical$gene),
         legend = "top", top = 0,
         select.top.method = "padj",
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = 12))

summary(temp34.regiontropical)
head(temp34.regiontropical[order(temp34.regiontropical$padj), ])

#Tropical vs. Temperate Interaction @ 36°C
ggmaplot(temp36.regiontropical, main = "Tropical vs. Temperate @ 36°C",
         fdr = 0.05, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(temp36.regiontropical$gene),
         legend = "top", top = 0,
         select.top.method = "padj",
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = 12))

summary(temp36.regiontropical)
head(temp36.regiontropical[order(temp36.regiontropical$padj), ])

#PCA
vsd <- vst(ddslrtemppop, blind=FALSE)
data <- plotPCA(vsd, ntop= 29868, intgroup=c("temp","region"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=temp, shape=region)) +
  geom_point(size=4, alpha=0.85) +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values = c("#FEE0D2", "#FC9272", "#EF3B2C", "#99000D")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic() +
  theme(aspect.ratio=1)


#Venn Diagram
#venn diagram
temp_32_vs_25V<- rownames(subset(temp_32_vs_25, padj < 0.05))
temp_34_vs_25V<- rownames(subset(temp_34_vs_25, padj < 0.05))
temp_36_vs_25V<- rownames(subset(temp_36_vs_25, padj < 0.05))
region_tropical_vs_temperateV<- rownames(subset(region_tropical_vs_temperate, padj < 0.05))
temp32.regiontropicalV<- rownames(subset(temp32.regiontropical, padj < 0.05))
temp34.regiontropicalV<- rownames(subset(temp34.regiontropical, padj < 0.05))
temp36.regiontropicalV<- rownames(subset(temp36.regiontropical, padj < 0.05))

x <- rownames(subset(temp_32_vs_25, padj < 0.05))
y <- rownames(subset(temp_34_vs_25, padj < 0.05))
z <- rownames(subset(temp_36_vs_25, padj < 0.05))
a <- rownames(subset(region_tropical_vs_temperate, padj < 0.05))
b <- rownames(subset(temp32.regiontropical, padj < 0.05))
c <- rownames(subset(temp34.regiontropical, padj < 0.05))
d <- rownames(subset(temp36.regiontropical, padj < 0.05))

library(VennDiagram)
VennDiagram::venn.diagram(x = list(T32 = temp_32_vs_25V, 
                                   T34 = temp_34_vs_25V, 
                                   T36 = temp_36_vs_25V,
                                   Troptemp = region_tropical_vs_temperateV
                                   ),
                          category.names = c("T32°C", 
                                             "T34°C",
                                             "T36°C",
                                             "Regions"
                                             ),
                          imagetype = "png",
                          filename = "venn_diagram_interactions.png",
                          fill = c("blue", "red", "purple", "green"), 
                          alpha = 0.5,
                          label.col = "white", 
                          fontfamily = "serif", 
                          fontface = "bold",
                          cat.cex = 2,
                          cat.dist = .1,
                          margin = .1,
                          cex = 2,
                          cat.col = c("blue", "red", "purple", "green"), 
                          cat.fontfamily = "serif",
                          cat.fontface = "bold",
                          na = "stop")

#Heat Maps


#Testers 
# annotation_colors = list(
#   region = c(temperate="#7CAE00", tropical="#00BFC4"))


d <-plotCounts(dds, gene="FBtr0080714", intgroup = (c("pop","temp", "region")), returnData=TRUE)

#dlog <- cbind(d,log2(d$count))

#write.csv(d, "/Users/lockwoodlab/Desktop/RNA Seq//REG_sesB_RB.csv")

p <-ggplot(d, aes(x=temp, y=count, shape=region, colour = pop)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("25","32","34", "36")) +
  scale_shape_manual(name = "region", values = c(1, 2)) 
p + theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15)) +
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=15)) +
  labs(x="Temperature (°C)", y = "Normalized Transcript Abundance", shape = "Region", color = "Population", shape = "Region") #why shape not work?

#boxplot
p<-ggplot(d, aes(x=temp, y=count, fill=region)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#00BC59", "#00B8E5"),
                    label = c("Temperate", "Tropical"),
                    name = "Region") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15)) +
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=15)) +
  labs(x="Temperature (°C)", y = "Normalized Transcript Abundance")

p

#heat map factor order
library("pheatmap")
#top 100 transcripts 
vsd <- vst(ddslrttempreg, blind = FALSE, nsub = 10000)
vsd <- vst(ddslrttemp, blind = FALSE, nsub = 10000)
vsd <- vst(ddslrtreg, blind = FALSE, nsub = 10000)
#topgenes <- head(rownames(region_tropical_vs_temperate), 100)
#for region
newdata <- region_tropical_vs_temperate[ which(region_tropical_vs_temperate$padj < 0.05), ]
toppvalgenes <- rownames(newdata)

toppvalgenes <- rownames(region_tropical_vs_temperate)[order(region_tropical_vs_temperate, region_tropical_vs_temperate$padj, decreasing = TRUE)][1:798] #full

toppvalgenes <- rownames(region_tropical_vs_temperate)[order(region_tropical_vs_temperate, region_tropical_vs_temperate$padj, decreasing = TRUE)][1:828] #reduced

# #for temperature
toppvalgenes <- rownames(temp_32_vs_25)[order(temp_32_vs_25, temp_32_vs_25$padj, decreasing = TRUE)][1:4136]

# just TEMP
newdata <- temp_32_vs_25[ which(temp_32_vs_25$padj < 0.05), ]
toppvalgenes <- rownames(newdata)
#toppvalgenesREG <- rownames(newdata)

#toppvalgenes <- rownames(temp_32_vs_25)[order(temp_32_vs_25, temp_32_vs_25$padj, decreasing = FALSE)][1:4534]

#toppvalgenes <- rownames(temp_32_vs_25)[order(temp_32_vs_25, temp_32_vs_25$padj)][1:4534]

# toppvalgenes <- rownames(res_pop_CP_vs_BO)[order(res_pop_CP_vs_BO, res_pop_CP_vs_BO$padj, decreasing = TRUE)][1:100]
mat <- assay(vsd)[toppvalgenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(ddslrtreg)[, c("temp", "region")])
library("RColorBrewer")
#colors<-colorRampPalette(rev(brewer.pal(n=7,name="RdYlGn")))(255)
#my_palette <- colorRampPalette(c("dodgerblue", "black", "yellow"))(n = 300)
#my_palette <- c("green",colorRampPalette(colors = c("red", "black", "green")), n = 3)
# pheatmap(mat, 
#          annotation_col = df, 
#          show_colnames = FALSE, show_rownames = FALSE, cluster_rows=TRUE, cluster_cols=FALSE, fontsize_row = 7, col=my_palette, breaks=breaksList, border_color=NA)

#sets scale bins
breaksList = seq(-1, 1, by = 0.005)

#sets colors for factors
annotation_colors = list(
 region = c(temperate="springgreen3", tropical="cyan2"),
temp = c(`25`="#FEE0D2", `32`="#FC9272",`34`="#EF3B2C", `36`="#99000D"))
  
pheatmap(mat, annotation_col = df, annotation_colors = annotation_colors,
         show_colnames = FALSE, show_rownames = FALSE, cluster_rows=TRUE, cluster_cols=FALSE, clustering_distance_rows = "correlation", fontsize_row = 7, col=colorRampPalette(c("dodgerblue", "black", "yellow"))(length(breaksList)), breaks=breaksList, border_color=NA)

#



#brewer.pal(n = 8, name = "Reds")

 #changing transcript IDs to Gene IDs

x <- read_delim("dmel-all-r6.34.gtf", 
                delim = "\t", 
                col_names = c("chr", "source", "feature", "start", 
                              "end", "score", "strand", "frame", "attributes"))

x$attributes <- str_replace_all(x$attributes, "\"", "")

x <- x %>%
  separate(attributes, 
           sep = ";", 
           into = c("gene_id", "gene_symbol", "transcript_id", 
                    "transcript_symbol", "notes", "extra", "extraextra"))
x$gene_id <- str_replace_all(x$gene_id, "gene_id ", "")
x$gene_symbol <- str_replace_all(x$gene_symbol, " gene_symbol ", "")
x$transcript_id <- str_replace_all(x$transcript_id, " transcript_id ", "")
x$transcript_symbol <- str_replace_all(x$transcript_symbol, " transcript_symbol ", "")

# transcript_to_gene <- function (trans_vec, gtf_dat){
#   gtf_dat <- gtf_dat %>%
#     distinct(transcript_id, .keep_all = TRUE) %>%
#     filter(transcript_id != "")
#   dat <- tibble::enframe(trans_vec, name = NULL, value = "transcript_id")
#   dat$gene <- dat$transcript_id
#   for (i in 1:nrow(dat)){
#     print(i)
#     temp <- which(dat$transcript_id[i] == gtf_dat$transcript_id, TRUE)
#     if (length(temp) == 0) { 
#         gene_replace <- NA
#       } else {
#         gene_replace <- gtf_dat$gene_symbol[temp]
#       }
#     dat$gene <- gsub(dat$transcript_id[i], 
#                      gene_replace, 
#                      dat$gene)
#   }
#   
#   dat_wo_na <- dat$gene[-which(is.na(dat$gene))]
#   return(dat_wo_na)
# }


transcript_to_gene <- function (trans_vec, gtf_dat){
  gtf_dat <- gtf_dat %>%
    distinct(transcript_id, .keep_all = TRUE) %>%
    filter(transcript_id != "")
  dat <- tibble::enframe(trans_vec, name = NULL, value = "transcript_id")
  dat$gene <- dat$transcript_id
  for (i in 1:nrow(dat)){
    print(i)
    temp <- which(dat$transcript_id[i] == gtf_dat$transcript_id, TRUE)
    if (length(temp) == 0) {
      gene_replace <- NA
    } else {
      gene_replace <- gtf_dat$gene_symbol[temp]
    }
    dat$gene <- gsub(dat$transcript_id[i],
                     gene_replace,
                     dat$gene)
  }
  dat_wo_na <- dat$gene[which(!is.na(dat$gene))]
  return(dat_wo_na)
}


#rntemp_32_vs_25 <-read.csv("temp_32_vs_25_tempreg.csv")
#rntemp_36_vs_25 <-read.csv("temp_36_vs_25_tempreggenes.csv")


#gntempreg <- transcript_to_gene(rntemp_32_vs_25$X, x)
#gntempreg36temp <- transcript_to_gene(temp_36_vs_25_genes, x)
gntempreg_reg <- transcript_to_gene(region_tropical_vs_temperate_genes, x)



gntempreg <- write.csv(unique(gntempreg), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegntempreg.csv")
gntempreg36temp  <- write.csv(unique(gntempreg36temp), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegntempreg36temp.csv")
gntempreg_reg  <- write.csv(unique(gntempreg_reg), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegntempreg_reg.csv")

#temp/pop interaction
gntemppop_int <- transcript_to_gene(res_temp34.popVT9_genes, x)
gntemppop_int <- write.csv(unique(gntemppop_int), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegntemppop_int.csv")

#Temp + pop
gntemppop_pop <- transcript_to_gene(res_pop_CP_vs_BO_genes, x)
gntemppop_pop <- write.csv(unique(gntemppop_pop), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegntemppop_pop.csv")


#TEMP
gnTEMP <- transcript_to_gene(TEMP_temp_32_vs_25_genes, x)
gnTEMP <- write.csv(unique(gnTEMP), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegnTEMP.csv")

#1/8 options
#upall
gnUPall <- transcript_to_gene(Upall_temp, x)
write.csv(unique(gnUPall), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegnUpall.csv")

#downall
gndownall <- transcript_to_gene(downall_temp, x)
write.csv(unique(gndownall), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegndownall.csv")

#up 32 down rest
gnup32 <- transcript_to_gene(up_32_down_rest_temp, x)
write.csv(unique(gnup32), "/Users/lockwoodlab/Desktop/RNA Seq//gnup32.csv")

#down32 up rest 
gndown32 <- transcript_to_gene(down32, x)
write.csv(unique(gndown32), "/Users/lockwoodlab/Desktop/RNA Seq/gndown32.csv")

#up 34 down rest
gnup34 <- transcript_to_gene(up34, x)
write.csv(unique(gnup34), "/Users/lockwoodlab/Desktop/RNA Seq/gnup34.csv")

#up 36 down rest
gnup36 <- transcript_to_gene(up36, x)
write.csv(unique(gnup36), "/Users/lockwoodlab/Desktop/RNA Seq/gnup36.csv")

#down 36 up rest
gndown36 <- transcript_to_gene(down36, x)
write.csv(unique(gndown36), "/Users/lockwoodlab/Desktop/RNA Seq/gndown36.csv")

#REGION
gnREG <- transcript_to_gene(Reg_temp_32_vs_25_genes, x)
gnREG <- write.csv(unique(gnREG), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegnREG.csv")

gnREG <- transcript_to_gene(toppvalgenesREG, x)

#just trop
gnTROP <- transcript_to_gene(tropsigT, x)
gnTROP <- write.csv(unique(gnTROP), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegnTROP.csv")

gnTemperate <- transcript_to_gene(tempsigT, x)
gnTemperate <- write.csv(unique(gnTemperate), "/Users/lockwoodlab/Desktop/RNA Seq//uniquegnTemperate.csv")







