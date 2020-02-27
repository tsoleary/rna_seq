# tso_analysis.Rmd manahttan plot and other lost information -------------------


# require(AnnotationDbi)
# require(org.Dm.eg.db)
# require(ggcorrplot)


# # load all the gwas snps for the manhattan plot
# gwas_cold_manhat <- read_delim("CTmin/gwas.all.assoc", 
#                                delim = " ", 
#                                col_names = TRUE)
# gwas_hot_manhat <- read_delim("CTmax/gwas.all.assoc", 
#                               delim = " ", 
#                               col_names = TRUE)
# 
# # wrangle the manhat dfs for qqman
# gwas_cold_manhat <- gwas_cold_manhat %>%
#   separate(ID, into = c("CHR", "BP", "type")) %>%
#   mutate(CHR = case_when(CHR == "2L" ~ 1,
#                          CHR == "2R" ~ 2,
#                          CHR == "3L" ~ 3,
#                          CHR == "3R" ~ 4,
#                          CHR == "4" ~ 5,
#                          CHR == "X" ~ 6)) %>%
#   dplyr::filter(!is.na(AvgMixedPval) | !is.na(CHR) | !is.na(BP))
# 
# gwas_cold_manhat$P <- gwas_cold_manhat$AvgMixedPval
# gwas_cold_manhat$BP <- as.numeric(gwas_cold_manhat$BP)
# 
# 
# gwas_hot_manhat <- gwas_hot_manhat %>%
#   separate(ID, into = c("CHR", "BP", "type")) %>%
#   mutate(CHR = case_when(CHR == "2L" ~ 1,
#                          CHR == "2R" ~ 2,
#                          CHR == "3L" ~ 3,
#                          CHR == "3R" ~ 4,
#                          CHR == "4" ~ 5,
#                          CHR == "X" ~ 6)) %>%
#   dplyr::filter(!is.na(AvgMixedPval) | !is.na(CHR) | !is.na(BP))
# 
# gwas_hot_manhat$P <- gwas_hot_manhat$AvgMixedPval
# gwas_hot_manhat$BP <- as.numeric(gwas_hot_manhat$BP)


# setwd(here::here("DRGP_GWAS"))
# 
# png("gwas_cold_manhattan.png", width = 600, height = 300)
# qqman::manhattan(gwas_cold_manhat,
#                  annotateTop = FALSE,
#                  chrlabs = c("2L", "2R", "3L", "3R", "4", "X"),
#                  ylim = c(0, 8),
#                  suggestiveline = -log10(1e-04),
#                  genomewideline = -log10(1e-05))
# dev.off()
# 
# 
# png("gwas_hot_manhattan.png", width = 600, height = 300)
# qqman::manhattan(gwas_hot_manhat,
#                  annotateTop = FALSE,
#                  chrlabs = c("2L", "2R", "3L", "3R", "4", "X"),
#                  ylim = c(0, 8),
#                  suggestiveline = -log10(1e-04),
#                  genomewideline = -log10(1e-05))
# dev.off()


## Labeled by padj

**Cold v Ctrl MA plot**
  
  ```{r}
ggmaplot(res_cold, main = expression("4°C" %->% "25°C"),
         fdr = 0.01, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(res_cold$gene),
         legend = "top", top = 20,
         select.top.method = "padj",
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = legend.font.size)) 
```

**Hot v Ctrl MA plot**
  
  ```{r}
ggmaplot(res_hot, main = expression("37°C" %->% "25°C"),
         fdr = 0.01, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         genenames = as.vector(res_hot$gene),
         legend = "top", top = 20,
         select.top.method = "padj",
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = legend.font.size)) 
```


## GWAS labels

**Cold v Ctrl MA plot**
  
  ```{r}
deg_gwas_cold_distinct <- deg_gwas_cold %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::filter(baseMean > 0)

ggmaplot_test(deg_gwas_cold_distinct, main = expression("4°C" %->% "25°C"),
              fdr = 0.01, fc = 0, size = 1,
              color_sig = FALSE,
              genenames = as.vector(deg_gwas_cold_distinct$gene),
              legend = "top", top = 20,
              select.top.method = "gwas",
              font.label = c("bold", 11), 
              font.label.sig = c(11, "bold", "black"),
              label.rectangle = TRUE,
              font.legend = "bold",
              font.main = "bold",
              ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = legend.font.size)) + 
  guides(color = guide_legend(override.aes = list(shape = 20)))
```

**Hot v Ctrl MA plot**
  
  ```{r}
deg_gwas_hot_distinct <- deg_gwas_hot %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::filter(baseMean > 0)

ggmaplot_test(deg_gwas_hot_distinct, main = expression("37°C" %->% "25°C"),
              fdr = 0.01, fc = 0, size = 1,
              color_sig = FALSE,
              genenames = as.vector(deg_gwas_hot_distinct$gene),
              legend = "top", top = 11,
              select.top.method = "gwas",
              font.label = c("bold", 11), 
              font.label.sig = c(11, "bold", "black"),
              label.rectangle = TRUE,
              font.legend = "bold",
              font.main = "bold",
              ggtheme = ggplot2::theme_minimal()) + 
  theme(legend.text = element_text(size = legend.font.size)) + 
  guides(color = guide_legend(override.aes = list(shape = 20))) 
```


## GWAS scatter & density plots (10<sup>-4</sup>)

```{r, fig.height = 7, fig.width = 7}

plot_scatter_side_density_xy(deg_simple_distinct_4, 
                             x_ = "log2FoldChange.cold", 
                             y_ = "log2FoldChange.hot",
                             labs_x = "Cold vs Ctrl (Log2 fold-change)",
                             labs_y = "Hot vs Ctrl (Log2 fold-change)",
                             bg.string = "NS: 6485",
                             bg.density.string = "Other",
                             labs_sets_density = "Top GWAS-associated genes",
                             id_ = "id", 
                             set_ = "groupn", 
                             set_density_ = "gwas_g", 
                             sets.density.colors = c("red", "blue", "yellow"),
                             sets.colors = color_set,
                             n_auto_label = 0)
```

```{r, fig.height = 4, fig.width = 7}
ggplot(deg_simple_distinct_4) +
  geom_line(mapping = aes(x = log2FoldChange.hot, color = gwas_g), 
            data = deg_simple_distinct_4, stat = "density") +
  geom_line(mapping = aes(x = log2FoldChange.cold, color = gwas_g), 
            linetype = "dashed", data = deg_simple_distinct_4, stat = "density") +
  scale_color_manual(values = c("red", "blue", "grey50"), drop = FALSE) +
  labs(x = "log 2 fold change vs. Ctrl", 
       color = "GWAS", 
       caption = "Cold dashed line\nHot solid line") +
  theme_classic()
```