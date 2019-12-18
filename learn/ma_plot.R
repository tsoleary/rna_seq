# MA plot of RNAseq data for entire dataset ------------------------------------
# https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html

library(ggpubr)

directory_results <- "/Users/tsoleary/R/rna_seq/results"
directory_plots <- "/Users/tsoleary/R/rna_seq/plots"
  
# groups being compared as how it was saved in the result folder
compCond <- "hot_cold"

# set outputPrefix
prefix <- "Dm_DESeq2"
outputPrefix <- paste(prefix, compCond, sep = "_")

# load result file
setwd(directory_results)
csv_result_file <- list.files()[grepl(compCond, list.files())]
res <- read.csv(csv_result_file)

# create ma plot
ggmaplot(res, main = compCond,
         fdr = 0.05, fc = 2, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(res$gene),
         legend = "top", top = 20,
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())


setwd(directory_plots)
dev.copy(png, paste0(outputPrefix, "_MAplot_ggpubr_analysis.png"))
dev.off()
