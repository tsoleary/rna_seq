# Cahan Lab DEG ----------------------------------------------------------------

require(tidyverse)
require(DESeq2)
require(ggpubr)
require(VennDiagram)
require(AnnotationDbi)
require(org.Dm.eg.db)

directory_counts <- here::here("cahan/counts")
directory_results <- here::here("cahan/results")
directory_plots <- here::here("cahan/plots")
directory_gwas <- here::here("DRGP_GWAS")

prefix <- "Dm_cahan_deg"

# Import in .counts files ------------------------------------------------------

# working directory containg the .counts files and metadata file
setwd(directory_counts)
metadata <- read.csv("chaos.txt")

# read in each .cnt and compile a dataframe of read counts for each sample
count_df <- data.frame()

for (i in 1:length(as.character(metadata$file))){
  t <- read.table(as.character(metadata$file)[i], header = TRUE)
  t <- t[, c(1, length(t))] # use only the first and last column (gene & counts)
  colnames(t) <- c("gene", as.character(metadata$label)[i])
  if (i == 1){
    count_df <- t 
  } else{
    count_df <- full_join(count_df, t, by = "gene")
  }
}

# create .csv file to upload to PIVOT for DESeq normalization
write.csv(count_df, paste0(prefix, "_counts.csv"), row.names = FALSE)

# convert to a matrix and save for DESeqDataSetFromMatrix
counts_mat <- as.matrix(count_df[,2:length(count_df)])
rownames(counts_mat) <- count_df$gene

# DESeq set up -----------------------------------------------------------------

treatments <- as.character(unique(metadata$condition))

ddsHTSeq <- DESeqDataSetFromMatrix(countData = counts_mat,
                                   colData = metadata,
                                   design = ~ condition)

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

# DESeq2 normalization and differential gene expression analysis ---------------
setwd(directory_results)

# normalization
dds <- DESeq(ddsHTSeq)

# set the ctrl as the group to be compared against
dds$condition <- relevel(dds$condition, ref = "ctrl")

# write normalized results 
write.table(as.data.frame(counts(dds, normalized = TRUE)), 
            file = paste0(prefix, "_norm_counts.txt"), sep = '\t')

# get the results of the differential expression
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
write.csv(resdata_hot, 
          file = paste0(prefix, "_hot", "_results_with_norm.csv"), 
          row.names = FALSE)
write.csv(resdata_cold, 
          file = paste0(prefix, "_cold", "_results_with_norm.csv"), 
          row.names = FALSE)

# MA plot of RNAseq data for entire dataset ------------------------------------

# groups being compared as how it was saved in the result folder
comp_cond <- "cold"
output_prefix <- paste(prefix, comp_cond, sep = "_")

# load result file
setwd(directory_results)
csv_result_file <- list.files()[grepl(comp_cond, list.files())]
res <- read.csv(csv_result_file)

# create ma plot
ggmaplot(res, main = expression("4°C" %->% "25°C"),
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
dev.copy(png, paste0(output_prefix, "_MAplot_ggpubr_analysis.png"))
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

# make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main = "Scree Plot", 
        xlab = "Principal Component", 
        ylab = "Percent Variation")

# now make a fancy looking plot that shows the PCs and the variation:
pca.data <- data.frame(Sample = rownames(pca$x),
                       X = pca$x[,1],
                       Y = pca$x[,2])

pca.data$group <- gsub("_[[:alnum:]]", "", pca.data$Sample)

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

# all differentially expressed genes
hot_deg <- row.names(subset(res_hot, res_hot$padj <= 0.05))
cold_deg <- row.names(subset(res_cold, res_cold$padj <= 0.05))

setwd(directory_plots)
venn.diagram(x = list(HOT = hot_deg, COLD = cold_deg),
             filename = paste0(prefix, "_venn_diagram.tiff"),
             lwd = 4, 
             fill = c("red", "blue"), 
             alpha = 0.5,
             label.col = "black", 
             cex = 1.7,
             fontfamily = "serif", 
             fontface = "bold",
             cat.col = c("red", "blue"), 
             cat.cex = 3, 
             cat.fontfamily = "serif",
             cat.fontface = "bold", 
             cat.dist = c(0.05, 0.03), 
             cat.pos = c(220, 160))

# subset to only genes that are significantly upregulated, L2FC > 1
hot_deg_up <- row.names(subset(res_hot, res_hot$padj <= 0.05 & 
                               res_hot$log2FoldChange > 1))
cold_deg_up <- row.names(subset(res_cold, res_cold$padj <= 0.05 &
                                res_cold$log2FoldChange > 1))

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

# subset to only genes that are significantly downregulated, L2FC < -1.5
hot_deg_down <- row.names(subset(res_hot, res_hot$padj <= 0.05 & 
                                 res_hot$log2FoldChange < -1))
cold_deg_down <- row.names(subset(res_cold, res_cold$padj <= 0.05 &
                                  res_cold$log2FoldChange < -1))

setwd(directory_plots)
venn.diagram(x = list(HOT = hot_deg_down, COLD = cold_deg_down),
             filename = paste0(prefix, "_down_reg_venn_diagram.tiff"),
             lwd = 4, 
             fill = c("red", "blue"), 
             alpha = 0.5,
             label.col = "black", 
             cex = 1.5,
             fontfamily = "serif", 
             fontface = "bold",
             cat.col = c("red", "blue"), 
             cat.cex = 3, 
             cat.fontfamily = "serif",
             cat.fontface = "bold", 
             cat.dist = c(0.05, 0.03), 
             cat.pos = c(220, 160))

# GWAS DGRP --------------------------------------------------------------------

# cold =======
outputPrefix <- paste0(prefix, "_cold")
# import GWAS data
setwd(directory_gwas)
gwas <- read.table("CTmin/gwas.all.assoc", header = TRUE)

# import DEG data
setwd(directory_results)
deg <- read.csv(list.files()[grepl("cold", list.files())])

# get the FBgn# for all the snps
gwas$FBgn <- str_extract(gwas$GeneAnnotation, 
                         "FBgn[[:digit:]]+") 

# FBgn to gene_symbol
gwas$gene <- as.character(mapIds(org.Dm.eg.db, 
                                 keys = gwas$FBgn, 
                                 column = "SYMBOL", 
                                 keytype = "FLYBASE",
                                 multiVals = "first")) 
comb <- full_join(deg, gwas, by = "gene")

# get the median or average for each treatment 
comb_avg <- comb %>% 
  dplyr::select(contains("_"), gene, ID, padj, log2FoldChange) %>%
  pivot_longer(contains("_"), names_to = "group", values_to = "expression") %>%
  mutate(group = str_replace(group, "_[[:digit:]]*$", "")) %>%
  group_by(gene, group, ID, padj, log2FoldChange) %>%
  summarize(expression = median(expression, na.rm = TRUE)) %>%
  filter(expression > 0) %>%
  pivot_wider(names_from = group, values_from = expression)

comb_sort <- comb_avg %>%
  mutate(g = case_when(is.na(ID) & padj < 0.05 & abs(log2FoldChange) > 1 ~ "DEG",
                       is.na(ID) & padj > 0.05 & abs(log2FoldChange) > 1 ~ "all",
                       is.na(ID) & abs(log2FoldChange) < 1 ~ "all",
                       !is.na(ID) & padj < 0.05 & 
                         abs(log2FoldChange) > 1 ~ "GWAS-DEG",
                       !is.na(ID) & padj > 0.05 & 
                         abs(log2FoldChange) > 1 ~ "GWAS-NS",
                       !is.na(ID) & abs(log2FoldChange) < 1 ~ "all")) %>%
  arrange(g) 

ggplot(comb_sort, aes(x= ctrl, y = cold)) + 
  geom_point(aes(x = ctrl, y = cold, color = g), alpha = 0.7) +
  scale_color_manual(values = c("#999999", "#E69F00", "#ff0000", "#000000"), 
                     breaks = c("DEG", "GWAS-DEG", "GWAS-NS")) +
  scale_x_log10() +
  scale_y_log10() + 
  ylab("Cold (log10 expression)") +
  xlab("Control (log10 expression)") +
  theme_classic() +
  theme(legend.title = element_blank())

setwd(directory_plots)
dev.copy(png, paste0(outputPrefix, "_expression_and_gwas.png"))
dev.off()

# hot =======
outputPrefix <- paste0(prefix, "_hot")

# import GWAS data
setwd(directory_gwas)
gwas <- read.table("CTmax/gwas.top.annot", header = TRUE)

# import DEG data
setwd(directory_results)
deg <- read.csv(list.files()[grepl("hot", list.files())])

# get the FBgn# for all the snps
gwas$FBgn <- str_extract(gwas$GeneAnnotation, 
                         "FBgn[[:digit:]]+") 

# FBgn to gene_symbol
gwas$gene <- as.character(mapIds(org.Dm.eg.db, 
                                 keys = gwas$FBgn, 
                                 column = "SYMBOL", 
                                 keytype = "FLYBASE",
                                 multiVals = "first")) 
comb <- full_join(deg, gwas, by = "gene")

# get the median or average for each treatment 
comb_avg <- comb %>% 
  dplyr::select(contains("_"), gene, ID, padj, log2FoldChange) %>%
  pivot_longer(contains("_"), names_to = "group", values_to = "expression") %>%
  mutate(group = str_replace(group, "_[[:digit:]]*$", "")) %>%
  group_by(gene, group, ID, padj, log2FoldChange) %>%
  summarize(expression = median(expression, na.rm = TRUE)) %>%
  filter(expression > 0) %>%
  pivot_wider(names_from = group, values_from = expression)

# create discrete groups for each category
comb_sort <- comb_avg %>%
  mutate(g = case_when(is.na(ID) & padj < 0.05 & abs(log2FoldChange) > 1 ~ "DEG",
                       is.na(ID) & padj > 0.05 & abs(log2FoldChange) > 1 ~ "all",
                       is.na(ID) & abs(log2FoldChange) < 1 ~ "all",
                       !is.na(ID) & padj < 0.05 & 
                         abs(log2FoldChange) > 1 ~ "GWAS-DEG",
                       !is.na(ID) & padj > 0.05 & 
                         abs(log2FoldChange) > 1 ~ "GWAS-NS",
                       !is.na(ID) & abs(log2FoldChange) < 1 ~ "all")) %>%
  arrange(g) 

# count the number of genes in each category
count_df <- comb_sort %>%
  group_by(g) %>%
  summarize(count = n())

# plot the expression of the ctrl and hot/cold with each group labeled
ggplot(comb_sort, aes(x = ctrl, y = hot)) + 
  geom_point(aes(x = ctrl, y = cold, color = g), alpha = 0.7) +
  scale_color_manual(values = c("#999999", "#E69F00", "#ff0000", "#000000"), 
                     breaks = c("DEG", "GWAS-DEG", "GWAS-NS")) +
  scale_x_log10() +
  scale_y_log10() + 
  ylab("Hot (log10 expression)") +
  xlab("Control (log10 expression)") +
  theme_classic() +
  theme(legend.title = element_blank())

# save the plot
setwd(directory_plots)
dev.copy(png, paste0(outputPrefix, "_expression_and_gwas.png"))
dev.off()



# gwas brent stuff -------------------------------------------------------------

# seth unique plot

# import DEG data
setwd(directory_results)
deg_cold <- read.csv(list.files()[grepl("cold", list.files())])
deg_hot <- read.csv(list.files()[grepl("hot", list.files())])


deg <- full_join(deg_cold, deg_hot, by = "gene", suffix = c(".cold", ".hot"))

deg_grouped <- deg %>%
  mutate(cold_genes = case_when(padj.cold < 0.05 ~ "sig",
                                padj.cold >= 0.05 ~ "ns"),
         hot_genes = case_when(padj.hot < 0.05 ~ "sig",
                               padj.hot >= 0.05 ~ "ns"),
         cold_fc = case_when(log2FoldChange.cold < 0 ~ "down",
                             log2FoldChange.cold > 0 ~ "up"),
         hot_fc = case_when(log2FoldChange.hot < 0 ~ "down",
                            log2FoldChange.hot > 0 ~ "up")) %>%
  filter(!is.na(cold_genes) & !is.na(hot_genes) &
           !is.na(cold_fc) & !is.na(hot_fc)) %>%
  mutate(col = paste(cold_genes, hot_genes, cold_fc, hot_fc)) %>%
  mutate(group = case_when(col == "sig sig up up" |
                             col == "sig sig down down" ~ "shared",
                           col == "sig ns up up" |
                             col == "ns sig up up" |
                             col == "sig ns down down" |
                             col == "ns sig down down" ~ "ns",
                           col == "ns ns down down" |
                             col == "ns ns down up" |
                             col == "ns ns up down" |
                             col == "ns ns up up" ~ "ns",
                           col == "sig sig down up" |
                             col == "sig sig up down" ~ "unique",
                           col == "ns sig up down" |
                             col == "ns sig down up" |
                             col == "sig ns down up" |
                             col == "sig ns up down" ~ "ns")) %>%
  add_count(group) %>% 
  mutate(groupn = paste0(group, ' (', n, ')')) %>%
  arrange(groupn)
            

ggplot(deg_grouped, 
       aes(x = log2FoldChange.cold, y = log2FoldChange.hot, color = groupn)) +
  geom_point() + 
  scale_color_manual(values = c("#999999", "#E69F00", "#ff0000")) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  ylab("Hot vs Ctrl") +
  xlab("Cold vs Ctrl") + 
  theme_bw()

