# ------------------------------------------------------------------------------
# GRFP Pre-liminary data figures
# September 12, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Load packages
require(tidyverse)
require(DESeq2)
source(here::here("functions.R"))

# Load data
lt_df <- read_delim("~/Downloads/lockwood_2018_LT50_data.txt", 
                    delim = "\t")

# LT50 Phenotype data

# Figure A
p1 <- lt_df %>%
  filter(Region == "Tropical" | StateCountry == "Vermont_USA") %>%
  ggplot(aes(x = Region, y = Embryo_LT50, fill = Region)) +
  geom_boxplot() +
  geom_jitter(color = "black", fill = "grey50",
              shape = 21,
              size = 3,
              position = position_jitter(0.1)) +
  scale_fill_manual(values = c("dodgerblue", "firebrick3")) +
  ylab(expression("Embryo LT"[50]*" (°C)")) +
  labs(title = "Embryonic\nthermal tolerance") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.title.x = element_blank())

# Figure B
ddslrtreg <- readRDS(here::here("emily/counts/ddslrtreg.rds"))
Idh <- "FBtr0076670"
Ucp4A <- "FBtr0074464" 
Sod3_RA <- "FBtr0089938"
Sod3_RB <- "FBtr0089939"
sesB <- "FBtr0073421"
rl <- "FBtr0345337"


d <- plotCounts(ddslrtreg, 
                gene = Ucp4A, 
                intgroup = c("region","temp"), 
                returnData = TRUE)

p2 <- ggplot(d, aes(x = region, 
              y = count, 
              #fill = region
              fill = temp
              )) + 
  geom_boxplot() +
  # geom_jitter(color = "black", fill = "grey50",
  #             shape = 21,
  #             size = 1,
  #             position = position_jitter(0.1)) +
  scale_fill_manual(name = "Temperature",
                     labels = c("25°C", "32°C", "34°C", "36°C"),
                     values = RColorBrewer::brewer.pal(n = 4, name = "Paired")) +
  # scale_fill_manual(values = c("dodgerblue", "firebrick3")) +
  scale_x_discrete(limits = c("temperate", "tropical"),
                   labels = c("Temperate", "Tropical")) +
  labs(title = expression(paste(italic("Ucp4A"), " expression")), 
       y = "Normalized\ntrancript abundance") +
  theme_classic(base_size = 18) +
  #ylim(c(0, 800)) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  theme(axis.title.x = element_blank(),
        #legend.background = element_rect(linetype = 1, color = 1),
        legend.position = c(.3, .85),
        legend.key.height = unit(0.7, "line"),
        legend.key.width = unit(0.7, "line"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

cowplot::plot_grid(p1, p2, 
                   rel_widths = c(1.5, 2), 
                   labels = "AUTO")

ggsave(here::here("emily/scripts/grfp_fig_bs_18.png"), width = 8, height = 3.5)


# Por's informaiton
fst_ch <- read_csv(here::here("por/fst_CHxVT10_500_tidy_gene_annot.csv"))
fst_sk <- read_csv(here::here("por/fst_SKxVT8_500_tidy_gene_annot.csv"))

# Chiapas x VT10 cross

fst_ch_top1 <- fst_ch %>%
  filter(VT10vVT10F > quantile(VT10vVT10F, 0.9)) %>%
  filter(!is.na(gene_assoc)) 

num_genes <- max(stringr::str_count(fst_ch_top1$gene_assoc, pattern = ";") + 1)
y <- fst_ch_top1 %>%
  separate(gene_assoc, 
           into = paste("gene", 1:num_genes, sep = "_"), 
           sep = ";") %>%
  pivot_longer(contains("gene_"), 
               names_to = "lab", 
               values_to = "gene",
               values_drop_na = TRUE) %>% 
  distinct(gene) %>%
  mutate(g = "VT10F")

fst_ch_top2 <- fst_ch %>%
  filter(CHFvVT10 > quantile(CHFvVT10, 0.9)) %>%
  filter(!is.na(gene_assoc))

num_genes <- max(stringr::str_count(fst_ch_top2$gene_assoc, pattern = ";") + 1)
y2 <- fst_ch_top2 %>%
  separate(gene_assoc, 
           into = paste("gene", 1:num_genes, sep = "_"), 
           sep = ";") %>%
  pivot_longer(contains("gene_"), 
               names_to = "lab", 
               values_to = "gene",
               values_drop_na = TRUE) %>% 
  distinct(gene) %>%
  mutate(g = "CHF")

# St Kitts x VT8 cross
fst_sk_top1 <- fst_sk %>%
  filter(SKFvVT8 > quantile(SKFvVT8, 0.9)) %>%
  filter(!is.na(gene_assoc)) 

num_genes <- max(stringr::str_count(fst_sk_top1$gene_assoc, pattern = ";") + 1)
z <- fst_sk_top1 %>%
  separate(gene_assoc, 
           into = paste("gene", 1:num_genes, sep = "_"), 
           sep = ";") %>%
  pivot_longer(contains("gene_"), 
               names_to = "lab", 
               values_to = "gene",
               values_drop_na = TRUE) %>% 
  distinct(gene) %>%
  mutate(g = "SKF")

fst_sk_top2 <- fst_sk %>%
  filter(VT8vVT8F > quantile(VT8vVT8F, 0.9)) %>%
  filter(!is.na(gene_assoc))

num_genes <- max(stringr::str_count(fst_sk_top2$gene_assoc, pattern = ";") + 1)
z2 <- fst_sk_top2 %>%
  separate(gene_assoc, 
           into = paste("gene", 1:num_genes, sep = "_"), 
           sep = ";") %>%
  pivot_longer(contains("gene_"), 
               names_to = "lab", 
               values_to = "gene",
               values_drop_na = TRUE) %>% 
  distinct(gene) %>%
  mutate(g = "VT8F")

# Join all these together
deg_region <- read_csv(here::here("emily/results/tropvtempREGION.csv")) %>%
  filter(padj < 0.05)

deg_region <- transcript_to_gene_df(deg_region, gtf_df) %>%
  distinct(gene, .keep_all = TRUE)


gene_dat <- full_join(y, y2, by = "gene", suffix = c(".1", ".2")) %>%
  full_join(z, by = "gene", suffix = c("", ".3")) %>%
  full_join(z2, by = "gene", suffix = c(".3", ".4")) %>%
  mutate(g = paste(g.1, g.2, g.3, g.4, sep = "_")) %>%
  mutate(g = str_replace_all(g, "NA_", "")) %>%
  mutate(g = str_replace_all(g, "_NA", "")) %>%
  mutate(g_num = str_count(g, "_") + 1) %>%
  arrange(desc(g_num), g)

por_emily_df <- gene_dat %>%
  full_join(deg_region %>%
              mutate(g_exp = "RegionDEG") %>%
              select(gene, g_exp), by = "gene") %>%
  mutate(g = paste(g, g_exp, sep = "_")) %>%
  arrange(g_exp, desc(g_num), g) %>%
  mutate(g = str_replace_all(g, "NA_", "")) %>%
  mutate(g = str_replace_all(g, "_NA", "")) %>%
  mutate(g_num = str_count(g, "_") + 1) %>%
  arrange(g_exp, desc(g_num), g)

write_csv(por_emily_df, here::here("emily/results/exp_fst_gene_overlap.csv"))


por_emily_df %>%
  filter(gene %in% c("Atg6","Ccs","CG6762","dj-1beta",
                     "Gyc89Db","Idh","Keap1","Prx2540-1",
                     "rl","Sdhaf3","sesB","Sod3","Ucp4A")) %>%
  select(-c(g, g_num))

osr <- read_delim("~/Downloads/osr.txt", delim = " ")

test <- por_emily_df %>%
  filter(gene %in% osr$symbol) 
