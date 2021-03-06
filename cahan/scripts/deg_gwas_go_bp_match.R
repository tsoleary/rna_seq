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


# NEW ORA GO RESULTS with new gene set sizes -----------------------------------

setwd(here::here("cahan/results/GO_results_files"))

# load df and save columns and names we want
ora_CTmin <- read_tsv("CTmin_ORA_030120_enrichment_results.txt") %>%
  select(description, userId, FDR) %>%
  rename(GO = description, gene = userId)
ora_CTmax <- read_tsv("CTmax_ORA_030120_enrichment_results.txt") %>%
  select(description, userId, FDR) %>%
  rename(GO = description, gene = userId)

# load deg go results to see which match
ora_cold <- read_tsv("Cold_DEG_ORA_030120_enrichment_results.txt") 
ora_hot <- read_tsv("Hot_DEG_ORA_030120_enrichment_results.txt") 


# match the go cats to ones that appear in the deg go cats
ora_CTmin <- ora_CTmin %>%
  mutate(match = case_when(GO %in% ora_cold$description ~ "match",
                           !(GO %in% ora_cold$description) ~ "no")) %>%
  mutate(GO = paste(GO, match, FDR, sep = ";")) %>%
  select(GO, gene)

ora_CTmax <- ora_CTmax %>%
  mutate(match = case_when(GO %in% ora_hot$description ~ "match",
                           !(GO %in% ora_hot$description) ~ "no")) %>%
  mutate(GO = paste(GO, match, FDR, sep = ";")) %>%
  select(GO, gene)

# match with the pvals and lfc
ora_CTmin_fig <- match_p_lfc(dat = ora_CTmin, dat_ct = gwas_cold, dat_deg = res_cold)
ora_CTmax_fig <- match_p_lfc(dat = ora_CTmax, dat_ct = gwas_hot, dat_deg = res_hot)


# write df to txt file
setwd(here::here("cahan/results/GO_results_files"))
write_tsv(ora_CTmin_fig, "CTmin_go_bp_fig.txt")
write_tsv(ora_CTmax_fig, "CTmax_go_bp_fig.txt")


# table for supplement
setwd(here::here("cahan/results/GO_results_files"))

read_tsv("Cold_DEG_ORA_030120_enrichment_results.txt") %>%
  select(geneSet, description, FDR, enrichmentRatio, userId) %>%
  mutate(FDR = signif(FDR, 3)) %>%
  rename(`GO ID` = geneSet,
         `GO category` = description, 
         `Enrichment Ratio` = enrichmentRatio,
         Genes = userId) %>%
  mutate(Genes = str_replace_all(Genes, ";", ", ")) %>%
  write_tsv(here::here("cahan/results/GO_results_files/cold_deg_ora.txt"))

read_tsv("Hot_DEG_ORA_030120_enrichment_results.txt") %>%
  select(geneSet, description, FDR, enrichmentRatio, userId) %>%
  mutate(FDR = signif(FDR, 3)) %>%
  rename(`GO ID` = geneSet,
         `GO category` = description, 
         `Enrichment Ratio` = enrichmentRatio,
         Genes = userId) %>%
  mutate(Genes = str_replace_all(Genes, ";", ", ")) %>%
  write_tsv(here::here("cahan/results/GO_results_files/hot_deg_ora.txt"))


read_tsv("CTmin_ORA_030120_enrichment_results.txt") %>%
  select(geneSet, description, FDR, enrichmentRatio, userId) %>%
  mutate(FDR = signif(FDR, 3)) %>%
  rename(`GO ID` = geneSet,
         `GO category` = description, 
         `Enrichment Ratio` = enrichmentRatio,
         Genes = userId) %>%
  mutate(Genes = str_replace_all(Genes, ";", ", ")) %>%
  write_tsv(here::here("cahan/results/GO_results_files/ctmin_ora.txt"))

read_tsv("CTmax_ORA_030120_enrichment_results.txt") %>%
  select(geneSet, description, FDR, enrichmentRatio, userId) %>%
  mutate(FDR = signif(FDR, 3)) %>%
  rename(`GO ID` = geneSet,
         `GO category` = description, 
         `Enrichment Ratio` = enrichmentRatio,
         Genes = userId) %>%
  mutate(Genes = str_replace_all(Genes, ";", ", ")) %>%
  write_tsv(here::here("cahan/results/GO_results_files/ctmax_ora.txt"))



### write tsv for 



