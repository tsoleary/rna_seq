# gene set enrichment analysis -------------------------------------------------

#require(msigdbr)
require(clusterProfiler)
require(AnnotationDbi)
require(org.Dm.eg.db)
require(tidyverse)

# Set the prefix for all output file names
prefix <- "Dm_cahan_deg"

# Set the working directory
directory_results <- here::here("cahan/results")
setwd(directory_results)

deg_cold <- read.csv(list.files()[grepl("cold", list.files())])
deg_hot <- read.csv(list.files()[grepl("hot", list.files())])


deg <- deg_cold %>% 
  filter(deg_cold$padj < 0.05)

# how to convert symbol to another type supported by org.Dm.eg
gene.df <- bitr(as.character(deg$gene), 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Dm.eg.db,
                drop = FALSE)
head(gene.df)

deg <- full_join(deg, gene.df, by = c("gene" = "SYMBOL"))

# the groupGO defaults to the ENTREZID
# you can specify the type of key with keyType = "SYMBOL" but it wasn't working
# ggo <- groupGO(gene = as.character(gene.df$ENTREZID),
#                OrgDb = org.Dm.eg.db,
#                ont = "CC",
#                level = 3,
#                readable = TRUE)
# head(ggo)
# 
# ego2 <- enrichGO(gene = as.character(gene.df$ENTREZID),
#                  OrgDb = org.Dm.eg.db,
#                  ont = "CC",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff = 0.01,
#                  qvalueCutoff = 0.05,
#                  readable = TRUE)

# does the geneList have to be in ENTREZID format?





# the geneList has to have the a decending order of some vector


# gene list descending by significant padj
geneList <- deg %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::select(gene, padj) %>%
  dplyr::arrange(desc(padj))



# gene list descending by significant padj
geneList <- deg %>%
  dplyr::filter(padj < 0.05 & !is.na(ENTREZID)) %>%
  dplyr::select(ENTREZID, padj) %>%
  dplyr::arrange(desc(padj))

# gene list descending by abs(log_fc)
geneList <- deg %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::mutate(abs_lfc = abs(log2FoldChange)) %>%
  dplyr::select(gene, abs_lfc) %>%
  dplyr::arrange(desc(abs_lfc)) %>%
  dplyr::distinct(abs_lfc, .keep_all = TRUE)

geneListMFer <- as.numeric(geneList$abs_lfc) 
names(geneListMFer) <- geneList$gene


# ont	one of "BP", "MF", and "CC" subontologies, or "ALL" for all three
ego3 <- gseGO(geneList = geneListMFer,
              OrgDb = org.Dm.eg.db,
              ont = "CC",
              nPerm = 1000,
              minGSSize = 100,
              maxGSSize = 500,
              verbose = FALSE)
# i have no idea what the fuck is going on and imfuckingsickofit

# --> Expected input gene ID: 34974,31518,41836,40966,38327,36654
# Error in check_gene_id(geneList, geneSets) : 
#   --> No gene can be mapped....




