# frontiers paper --------------------------------------------------------------

# this script is to try to use the list that seth created that has the matched
# go bp categories for gwas and degs
# the goal is to have a list of go categories and the number of hits in each 
# category with the specific gene names in each line like it MacMillian et al 
# 2016 

# I am just going to focus on GO BP for CTmin & cold and CTmax & hot -----------
require(tidyverse)

# load the sheets with the true false values

setwd(here::here("cahan/results"))
bp_cold <- read_csv("COLD_GO_BP_CTmin.csv")
bp_hot <- read_csv("HOT_GO_BP_CTmax.csv")

x <- bp_cold %>%
  select_if(is.character) %>%
  dplyr::select(-ct_min) %>%
  pivot_longer(everything(),
               names_to = "cat",
               values_to = "gene") %>%
  filter(gene != "FALSE") %>%
  dplyr::select(cat:gene) %>%
  left_join(res_cold, by = "gene") %>%
  left_join(gwas_cold_all, by = "gene") %>%
  arrange(AvgMixedPval) %>%
  group_by(cat) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(cat, gene, log2FoldChange, padj, AvgMixedPval) %>%
  mutate(gene = paste(gene, log2FoldChange, padj, AvgMixedPval, sep = "_")) %>%
  dplyr::select(cat, gene) %>%
  pivot_wider(names_from = cat, 
              values_from = gene,
              values_fn = list(gene = list)) %>%
  pivot_longer(everything(),
               names_to = "cat",
               values_to = "gene_list")


# for loop it to unlist
x$genes <- vector(mode = "character", length = nrow(x))
for (i in 1:nrow(x)){
  x$genes[i] <- paste(unlist(x$gene_list[i]), sep = "", collapse = "_")
}

write_csv(x, "deg_gwas_cold_bp_table.csv")

# bp_hot
y <- bp_hot %>%
  select_if(is.character) %>%
  dplyr::select(-gene_CTmax) %>%
  pivot_longer(everything(),
               names_to = "cat",
               values_to = "gene") %>%
  filter(gene != "FALSE") %>%
  dplyr::select(cat:gene) %>%
  left_join(res_hot, by = "gene") %>%
  left_join(gwas_hot_all, by = "gene") %>%
  arrange(AvgMixedPval) %>%
  group_by(cat) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(cat, gene, log2FoldChange, padj, AvgMixedPval) %>%
  mutate(gene = paste(gene, log2FoldChange, padj, AvgMixedPval, sep = "_")) %>%
  dplyr::select(cat, gene) %>%
  pivot_wider(names_from = cat, 
              values_from = gene,
              values_fn = list(gene = list)) %>%
  pivot_longer(everything(),
               names_to = "cat",
               values_to = "gene_list")

# for loop it to unlist
y$genes <- vector(mode = "character", length = nrow(y))
for (i in 1:nrow(y)){
  y$genes[i] <- paste(unlist(y$gene_list[i]), sep = "", collapse = "_")
}


write_csv(y, "deg_gwas_hot_bp_table.csv")


# now lets do it but filter it out so it is only genes with AvgMixedPval < 10^-5 -----

x <- bp_cold %>%
  select_if(is.character) %>%
  dplyr::select(-ct_min) %>%
  pivot_longer(everything(),
               names_to = "cat",
               values_to = "gene") %>%
  filter(gene != "FALSE") %>%
  dplyr::select(cat:gene) %>%
  left_join(res_cold, by = "gene") %>%
  left_join(gwas_cold_all, by = "gene") %>%
  arrange(AvgMixedPval) %>%
  group_by(cat) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(cat, gene, log2FoldChange, AvgMixedPval) %>%
  dplyr::filter(AvgMixedPval < 10^-5) %>%
  mutate(gene = paste(gene, log2FoldChange, sep = "_")) %>%
  dplyr::select(cat, gene) %>%
  pivot_wider(names_from = cat, 
              values_from = gene,
              values_fn = list(gene = list)) %>%
  pivot_longer(everything(),
               names_to = "cat",
               values_to = "gene_list")


# for loop it to unlist
x$genes <- vector(mode = "character", length = nrow(x))
for (i in 1:nrow(x)){
  x$genes[i] <- paste(unlist(x$gene_list[i]), sep = "", collapse = "_")
}

write_csv(x, "deg_gwas_cold_bp_table-5.csv")

# bp_hot
y <- bp_hot %>%
  select_if(is.character) %>%
  dplyr::select(-gene_CTmax) %>%
  pivot_longer(everything(),
               names_to = "cat",
               values_to = "gene") %>%
  filter(gene != "FALSE") %>%
  dplyr::select(cat:gene) %>%
  left_join(res_hot, by = "gene") %>%
  left_join(gwas_hot_all, by = "gene") %>%
  arrange(AvgMixedPval) %>%
  group_by(cat) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(cat, gene, log2FoldChange, AvgMixedPval) %>%
  dplyr::filter(AvgMixedPval < 10^-5) %>%
  mutate(gene = paste(gene, log2FoldChange, sep = "_")) %>%
  dplyr::select(cat, gene) %>%
  pivot_wider(names_from = cat, 
              values_from = gene,
              values_fn = list(gene = list)) %>%
  pivot_longer(everything(),
               names_to = "cat",
               values_to = "gene_list")

# for loop it to unlist
y$genes <- vector(mode = "character", length = nrow(y))
for (i in 1:nrow(y)){
  y$genes[i] <- paste(unlist(y$gene_list[i]), sep = "", collapse = "_")
}


write_csv(y, "deg_gwas_hot_bp_table-5.csv")


# lets try joining the two
xy <- full_join(x, y, by = "cat")

write_csv(xy, "deg_gwas_hot_cold_comb_bp_table-5.csv")



# using the gene list that Brent has generated match up the go terms with genes 
# that are in the CTmin and CTmax ----------------------------------------------

setwd(here::here("cahan/results"))
ora_cold <- read_tsv("COLD_ORA_GOBP_genes.txt")
ora_hot<- read_tsv("HOT_ORA_GOBP_genes.txt")
gsea_cold <- read_tsv("COLD_GSEA_GOBP_genes.txt")
gsea_hot<- read_tsv("HOT_GSEA_GOBP_genes.txt")

ora_cold_match <- go_gene_match(dat = ora_cold, dat_match = gwas_cold)
ora_hot_match <- go_gene_match(dat = ora_hot, dat_match = gwas_hot)
gsea_cold_match <- go_gene_match(dat = gsea_cold, dat_match = gwas_cold)
gsea_hot_match <- go_gene_match(dat = gsea_hot, dat_match = gwas_hot)


write_tsv(ora_cold_match, "ora_cold_match_go_bp.txt")
write_tsv(ora_hot_match, "ora_hot_match_go_bp.txt")
write_tsv(gsea_cold_match, "gsea_cold_match_go_bp.txt")
write_tsv(gsea_hot_match, "gsea_hot_match_go_bp.txt")

