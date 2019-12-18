# venn diagram of DESeq results ------------------------------------------------

library(VennDiagram)

directory <- "/Users/tsoleary/R/rna_seq/plots"
setwd(directory)

# need dds object
# source("/Users/tsoleary/R/rna_seq/scripts/DESeq_r.R)

res_hot_con <- results(dds, contrast = c("sampleCondition", "HOT","CON"))
res_hot_con <- subset(res_hot_con, res_hot_con$padj <= 0.05)
hot_con_genes <- row.names(res_hot_con)

res_cold_con <- results(dds, contrast = c("sampleCondition", "COLD","CON"))
res_cold_con <- subset(res_cold_con, res_cold_con$padj <= 0.05)
cold_con_genes <- row.names(res_cold_con)

# Combining the two above..
comb <- c(hot_con_genes, cold_con_genes)

# Comparing comb with the above two
hot_con_genes_2 <- comb %in% hot_con_genes
cold_con_genes_2 <- comb %in% cold_con_genes 


# make a venn diagram
venn.diagram(x = list(HOT = hot_con_genes, COLD = cold_con_genes),
  filename = "deseq_venn_diagram.tiff",
  lwd = 4, fill = c("red", "blue"), alpha = 0.5,
  label.col = "white", cex = 3,
  fontfamily = "serif", fontface = "bold",
  cat.col = c("red", "blue"), cat.cex = 3, cat.fontfamily = "serif",
  cat.fontface = "bold", cat.dist = c(0.03, 0.03), cat.pos = c(200, 180))
