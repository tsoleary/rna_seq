# Pivot installation and dependencies ------------------------------------------
# TSO install October 10, 2019

# Install 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# dependecies that needs to be manually installed 

library("devtools")
BiocManager::install("GO.db")
BiocManager::install("HSMMSingleCell")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scater")

# Install PIVOT
install_github("qinzhu/PIVOT")
BiocManager::install("BiocGenerics") # You need the latest BiocGenerics >=0.23.3

library(PIVOT)
pivot()
