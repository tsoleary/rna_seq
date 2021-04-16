# ------------------------------------------------------------------------------
# Functions for Dan
# April 05, 2021
# TS O'Leary
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: exp_data_tidy
# Description: description
# Inputs: input_description
# Outputs: output_description

exp_data_tidy <- function(gene_fb, pheno_df) {
  exp %>%
    filter(gene == gene_fb) %>%
    pivot_longer(-gene, 
                 names_to = "sample", 
                 values_to = "norm_exp") %>%
    separate(sample, sep = ":", into = c("line", "rep", "sex")) %>%
    pivot_wider(names_from = sex, values_from = norm_exp) %>%
    pivot_wider(names_from = rep, values_from = c("male", "female")) %>%
    mutate(line = as.numeric(str_extract_all(line, "[0-9]+"))) %>%
    full_join(pheno_df, by = c("line" = "line_DGRP")) %>%
    pivot_longer(contains("ctm"), 
                 names_to = "type", 
                 values_to = "temp") %>%
    separate(type, sep = "_", into = c("limit", "sex_pheno")) %>%
    pivot_longer(contains("e_"), names_to = "type", values_to = "norm_exp") %>%
    separate(type, sep = "_", into = c("sex_exp", "rep")) %>%
    filter(sex_pheno == sex_exp)
} 
# End function -----------------------------------------------------------------



# ------------------------------------------------------------------------------
# Function: plot_mf_exp
# Description: Plots expression for males and females separately
# Inputs: input_description
# Outputs: output_description

plot_mf_exp <- function(dat) {
  ggplot(data = dat, 
         aes(x = temp, 
             y = norm_exp, 
             color = sex_exp)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", 
                se = FALSE) +
    labs(y = "Normalized Expression",
         x = "Critical Temperature (°C)") +
    scale_color_brewer(palette = "Dark2",
                       labels = c("Female", "Male"),
                       name = "Sex") +
    theme_classic()
} 
# End function -----------------------------------------------------------------



# ------------------------------------------------------------------------------
# Function: plot_min_max_gene
# Description: description
# Inputs: input_description
# Outputs: output_description

plot_min_max_gene <- function(gene_fb) {
  
  # Filter and tidy data for selected gene
  exp_g <- exp_data_tidy(gene_fb, pheno) %>%
    group_by(limit) %>%
    group_split()
  
  # Create plots for CT min and CTmax
  p_max <- plot_mf_exp(exp_g[[1]])
  p_min <- plot_mf_exp(exp_g[[2]]) + 
    theme(legend.position = "none")
  
  title <- cowplot::ggdraw() +
    cowplot::draw_label(paste(gene_fb, "expression"))
  
  # Plot both CTmin and CTmax
  plots <- cowplot::plot_grid(p_min, p_max, 
                              rel_widths = c(1, 1.25))
  
  cowplot::plot_grid(title, plots, nrow = 2, 
                     rel_heights = c(.1, 1))
  
} 
# End function -----------------------------------------------------------------