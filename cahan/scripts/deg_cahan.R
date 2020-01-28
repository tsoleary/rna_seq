# Cahan Lab DEG ----------------------------------------------------------------

require(tidyverse)
require(DESeq2)
require(ggpubr)
require(VennDiagram)
require(AnnotationDbi)
require(org.Dm.eg.db)
require(seqsetvis)
require(ssvRecipes)

source(here::here("functions.R"))

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



# seth unique plot -------------------------------------------------------------

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
                             col == "sig sig down down" ~ "Shared",
                           col == "sig ns up up" |
                             col == "ns sig up up" |
                             col == "sig ns down down" |
                             col == "ns sig down down" ~ "Sig1 Shared",
                           col == "ns ns down down" |
                             col == "ns ns down up" |
                             col == "ns ns up down" |
                             col == "ns ns up up" ~ "NS",
                           col == "sig sig down up" |
                             col == "sig sig up down" ~ "Unique",
                           col == "ns sig up down" |
                             col == "ns sig down up" |
                             col == "sig ns down up" |
                             col == "sig ns up down" ~ "Sig1 Unique")) %>%
  add_count(group) %>% 
  mutate(groupn = paste0(group, ': ', n)) %>%
  arrange(groupn)
            
clrs <- c("darkgrey", "goldenrod4", "goldenrod1", "coral", "coral4")

set.seed(1)
ggplot(deg_grouped, 
       aes(x = log2FoldChange.cold, y = log2FoldChange.hot, color = groupn)) +
  geom_point() + 
  scale_color_manual(values = clrs) +
  geom_hline(size = 1, yintercept = 0, color = "black") +
  geom_vline(size = 1, xintercept = 0, color = "black") +
  ylab("Hot vs Ctrl") +
  xlab("Cold vs Ctrl") + 
  theme_classic() + 
  labs(color = "") + 
  ggrepel::geom_label_repel(data = filter(deg_grouped, group == "Unique") %>%
                              top_n(20, wt = baseMean.cold),
                            aes(label = gene),
                            fill = "white",
                            color = 'black',
                            size = 3.5,
                            show.legend = FALSE) +
  theme(legend.position = "top")


# gwas brent stuff  ------------------------------------------------------------

# import GWAS data
setwd(directory_gwas)
gwas <- read.table("CTmax/gwas.top.annot", header = TRUE)
gwas_all <- read.table("CTmin/gwas.all.assoc", header = TRUE)

# import DEG data
setwd(directory_results)
deg_cold <- read.csv(list.files()[grepl("cold_results", list.files())])
deg_hot <- read.csv(list.files()[grepl("hot_results", list.files())])

deg <- full_join(deg_cold, deg_hot, by = "gene", suffix = c(".cold", ".hot"))


# Gene annotation 
# Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| 
# Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon [ | 
# ERRORS | WARNINGS ] )

# # convert fbgn to gene symbol -----
# # get the FBgn# for all the snps
# gwas$FBgn <- str_extract(gwas$GeneAnnotation,
#                          "FBgn[[:digit:]]+")
# # 
# # # FBgn to gene_symbol
# gwas$gene2 <- as.character(mapIds(org.Dm.eg.db,
#                                  keys = gwas$FBgn,
#                                  column = "SYMBOL",
#                                  keytype = "FLYBASE",
#                                  multiVals = "first"))
# # this is a worse method there is a way that some of the gene symbols 
# # are getting lost in translation

# split up the gene annotation line by the vertical pipe | ---------------------
annot_split <- str_split(gwas$GeneAnnotation, "\\|")

# the second item is the gene symbol
gwas$gene <- sapply(annot_split, "[[", 2)
# the third item is where the snp is located, intron, syn coding etc.
gwas$grps <- sapply(annot_split, "[[", 3)

# join the two data frames together
comb <- full_join(deg, gwas, by = "gene")

comb_filt <- comb %>%
  filter(AvgMixedPval < 0.00005)
 
ggplot(comb_filt, aes(x = abs(log(AvgMixedPval)), y = abs(log2FoldChange.hot))) +
  geom_point() 

ggplot(comb_filt, aes(x = AvgEff, y = -log(padj.hot))) +
  geom_point() + 
  ylim(0,25)

# create .rnk files for the GSEA in webgesault ---------------------------------

# import DEG data
setwd(directory_results)
deg_cold <- read.csv(list.files()[grepl("cold_results", list.files())])
deg_hot <- read.csv(list.files()[grepl("hot_results", list.files())])

# write files for gsea input
write.rnk(deg_cold, "cold_deg_p_01_gsea_lfc.rnk", padj_cut = 0.01)
write.rnk(deg_hot, "hot_deg_p_01_gsea_lfc.rnk", padj_cut = 0.01)



# import DEG data filter the unique and shared
# create a combined df with deg and gwas for each alone ------------------------
comb_cold <- full_join(res_cold, gwas_cold, by = "gene")
comb_hot <- full_join(res_hot, gwas_hot, by = "gene")


# create a data.frame with info for the individual ma_plots --------------------
deg_gwas_cold <- left_join(res_cold, gwas_cold, by = "gene")
deg_gwas_hot <- left_join(res_hot, gwas_hot, by = "gene")

# create a data.frame with all the info for the plotting later -----------------
temp <- left_join(deg, gwas_cold, by = "gene", suffix = c("", ".cold"))
total <- full_join(temp, gwas_hot, by = "gene", suffix = c(".cold", ".hot"))

# group by the significance
deg_grouped <- group_deg(total, pval_cut = 0.01, lfc_cut = 0)

uniq <- deg_grouped %>%
  dplyr::filter(group == "Unique")

uniq_sig1 <- deg_grouped %>%
  dplyr::filter(group == "Unique" | 
                  group == "Sig1 Unique")

shared_up <- deg_grouped %>%
  dplyr::filter(col == "sig sig up up" | 
                  col == "sig ns up up" |
                  col == "ns sig up up")
  


# write the .rnk files for gsea input
# or should I just do the over enrichment analysis

# write files for overenrichment analysis -----------
uniq %>%
  dplyr::select(gene) %>%
  write_delim("unique_p_01_gsea_lfc.txt", col_names = FALSE, delim = "\t")

uniq_sig1 %>%
  dplyr::select(gene) %>%
  write_delim("unique_sig_1_p_01_gsea_lfc.txt", col_names = FALSE, delim = "\t")

shared_up %>%
  dplyr::select(gene) %>%
  write_delim("shared_up_p_01_gsea_lfc.txt", col_names = FALSE, delim = "\t")


# logfoldchange expression scatter and density plot ----------------------------

# example call of custom plotting functions ------------------------------------
n = 50
xy_data = rbind(
  data.table(x = rnorm(10*n, 0, 1), y = rnorm(10*n, 0, 1), set = "background", set_density = "background"),
  data.table(x = rnorm(2*n, 2, 1), y = rnorm(2*n, 0, 1), set = "set1", set_density = "sets"),
  data.table(x = rnorm(2*n, 0, 1), y = rnorm(2*n, 2, 1), set = "set2", set_density = "sets"),
  data.table(x = rnorm(2*n, 2, 1), y = rnorm(2*n, 2, 1), set = "set3", set_density = "sets")
)
xy_data$id = seq_len(nrow(xy_data))

plot_scatter_side_density.xy(xy_data, x_ = "x", y_ = "y")
plot_scatter_side_density_xy(xy_data, x_ = "x", y_ = "y")
plot_scatter_side_density_xy(xy_data, x_ = "x", y_ = "y", 
                             set_density_ = "set_density",
                             sets.density.colors = c("red", "blue"))


# scatter and density plot of my data ------------------------------------------

color_set <- c("goldenrod4", "goldenrod1", "coral", "coral4", "grey50")

# the function from seqsetvis requires a sequential column in the data.frame
deg_grouped$id <- seq_len(nrow(deg_grouped))

deg_simple <- deg_grouped %>%
  dplyr::select(log2FoldChange.hot, log2FoldChange.cold, groupn, id, ID.hot, ID.cold) %>%
  mutate(gwas_g = case_when(is.na(ID.hot) & is.na(ID.cold) ~ "Other",
                            !is.na(ID.hot) & is.na(ID.cold)  ~ "CTmax",
                            is.na(ID.hot) & !is.na(ID.cold)  ~ "CTmin",
                            !is.na(ID.hot) & !is.na(ID.cold)  ~ "both"))


ggplot(deg_simple_distinct) +
  geom_line(mapping = aes(x = log2FoldChange.hot, color = gwas_g), 
            data = deg_simple_distinct, stat = "density") +
  geom_line(mapping = aes(x = log2FoldChange.cold, color = gwas_g), 
            linetype = "dashed", data = deg_simple_distinct, stat = "density") +
  scale_color_manual(values = c("red", "blue", "grey50"), drop = FALSE) +
  labs(x = "") +
  theme_classic()

p + ggplot() +
  geom_line( data = deg_simple, mapping = aes(x = log2FoldChange.cold, color = gwas_g), data = deg_simple, stat = "density")


# order the levels of the factor for plotting purposes to get NS on bottom
deg_simple$groupn <- as.factor(deg_simple$groupn)

deg_simple$groupn <- factor(deg_simple$groupn, levels = c("Shared: 6123", 
                                                          "Sig1 Shared: 3217", 
                                                          "Unique: 14",
                                                          "Sig1 Unique: 104",
                                                          "NS: 6437"))


plot_scatter_side_density_xy(test, 
                             x_ = "log2FoldChange.cold", 
                             y_ = "log2FoldChange.hot",
                             bg.string = "NS: 6437",
                             bg.density.string = "background",
                             labs_sets_density = "GWAS",
                             id_ = "id", 
                             set_ = "groupn", 
                             set_density_ = "gwas_g", 
                             sets.density.colors = c("yellow", "red", "blue"),
                             sets.colors = color_set,
                             n_auto_label = 0)


# feature mapper snps gwas ------

# load the gwas results
setwd(here::here("DRGP_GWAS"))
gwas_cold_all <- read.table("CTmin/gwas.top.annot", header = TRUE)
gwas_hot_all <- read.table("CTmax/gwas.top.annot", header = TRUE)

gwas_hot <- gwas_hot_all %>%
  dplyr::filter(AvgMixedPval < 10^-5)

gwas_cold <- gwas_cold_all %>%
  dplyr::filter(AvgMixedPval < 10^-5)

gwas_hot <- gwas_hot %>%
  separate(ID, c("chr", "pos", "type")) %>%
  dplyr::mutate(featureMap = paste(paste(chr, pos, sep = ":"), pos, sep = "-"))

write_delim(gwas_hot %>% dplyr::select(featureMap), 
            "gwas_hot_snps.txt", delim = "^t", col_names = FALSE)


gwas_cold <- gwas_cold %>%
  separate(ID, c("chr", "pos", "type")) %>%
  dplyr::mutate(featureMap = paste(paste(chr, pos, sep = ":"), pos, sep = "-"))

write_delim(gwas_cold %>% dplyr::select(featureMap), 
            "gwas_cold_snps.txt", delim = "^t", col_names = FALSE)

# snps for over-represenation analysis -----


setwd(here::here("cahan/results/ORA/ctmin"))
gwas_cold_all %>% 
  dplyr::select(gene) %>% 
  distinct(gene) %>% 
  filter(gene != "") %>% 
  write_delim("gwas_ctmin_genes_4_ora.txt", 
              delim = "^t", col_names = FALSE)


setwd(here::here("cahan/results/ORA/ctmax"))
gwas_hot_all %>% 
  dplyr::select(gene) %>% 
  distinct(gene) %>% 
  filter(gene != "") %>% 
  write_delim("gwas_ctmax_genes_4_ora.txt", 
              delim = "^t", col_names = FALSE)


# differentially expressed genes for over representation analysis --------------

# all differentially expressed genes
setwd(here::here("cahan/results/ORA/cold"))
res_cold %>% 
  dplyr::filter(padj < 0.01) %>%
  dplyr::select(gene) %>% 
  distinct(gene) %>% 
  filter(gene != "") %>% 
  write_delim("deg_cold_all_01_ora.txt", 
              delim = "^t", col_names = FALSE)

setwd(here::here("cahan/results/ORA/hot"))
res_hot %>% 
  dplyr::filter(padj < 0.01) %>%
  dplyr::select(gene) %>% 
  distinct(gene) %>% 
  filter(gene != "") %>% 
  write_delim("deg_hot_all_01_ora.txt", 
              delim = "^t", col_names = FALSE)


# only differentially expressed up
setwd(here::here("cahan/results/ORA/cold"))
res_cold %>% 
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  dplyr::select(gene) %>% 
  distinct(gene) %>% 
  filter(gene != "") %>% 
  write_delim("deg_cold_up_01_ora.txt", 
              delim = "^t", col_names = FALSE)

setwd(here::here("cahan/results/ORA/hot"))
res_hot %>% 
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  dplyr::select(gene) %>% 
  distinct(gene) %>% 
  filter(gene != "") %>% 
  write_delim("deg_hot_up_01_ora.txt", 
              delim = "^t", col_names = FALSE)

# only differentially expressed down
setwd(here::here("cahan/results/ORA/cold"))
res_cold %>% 
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  dplyr::select(gene) %>% 
  distinct(gene) %>% 
  filter(gene != "") %>% 
  write_delim("deg_cold_down_01_ora.txt", 
              delim = "^t", col_names = FALSE)

setwd(here::here("cahan/results/ORA/hot"))
res_hot %>% 
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  dplyr::select(gene) %>% 
  distinct(gene) %>% 
  filter(gene != "") %>% 
  write_delim("deg_hot_down_01_ora.txt", 
              delim = "^t", col_names = FALSE)






 
