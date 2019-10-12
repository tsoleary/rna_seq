# example gene count data ----

# https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf

library("DESeq2")
# BiocManager::install("parathyroidSE")
library("parathyroidSE")

data("parathyroidGenesSE")
se <- parathyroidGenesSE
colnames(se) <- se$run


library(PIVOT)
pivot()
