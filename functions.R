# functions for rna seq analysis -----------------------------------------------


# function to write .rnk files for gsea input ----------------------------------
write.rnk <- function(.data, file, padj_cut = 0.05){
  .data %>%
    filter(padj < padj_cut) %>%
    dplyr::select(gene, log2FoldChange) %>%
    filter(!is.na(log2FoldChange)) %>%
    arrange(log2FoldChange) %>%
    write_delim(file, col_names = FALSE, delim = "\t")
}


# custom plot_scatter_side_density_xy ------------------------------------------
plot_scatter_side_density_xy = function( xy_data,
                                         x_ = "x",
                                         y_ = "y",
                                         id_ = "id",
                                         set_ = "set",
                                         set_density_ = set_,
                                         #labels
                                         labs_x = x_,
                                         labs_y = y_,
                                         labs_sets = set_,
                                         labs_sets_density = set_density_,
                                         main_title = NULL,
                                         main_title.x = .02,
                                         main_title.y = .5,
                                         main_title.hjust = 0,
                                         main_title.vjust = .5,
                                         #point color and sizing
                                         sets.colors = NULL,
                                         sets.density.colors = NULL,
                                         bg.string = "background",
                                         bg.density.string = "background",
                                         bg.color = "gray70",
                                         bg.density.color = "gray70",
                                         sets.sizes = 1,
                                         sets.density.sizes = 1,
                                         bg.size = .5,
                                         bg.density.size = 0.5,
                                         #limits
                                         xlim_ = NULL,
                                         ylim_ = NULL,
                                         #point labelling
                                         n_auto_label = 8,
                                         manual_label = NULL,
                                         label_size = 6,
                                         label_color = 'black',
                                         label_GEOM = shadowtext::geom_shadowtext,
                                         #reference lines
                                         ref_line.x = 0,
                                         ref_line.x.color = "gray50",
                                         ref_line.y = 0,
                                         ref_line.y.color = "gray50",
                                         ref_line.slope = 1,
                                         ref_line.slope.color = "black",
                                         suppress_plot = FALSE){
  if(is.matrix(xy_data)){
    rn = rownames(xy_data)
    xy_data = as.data.table(xy_data)
    xy_data[[id_]] = rn
  }
  if(is.data.frame(xy_data) & !is.data.table(xy_data)){
    if(is.null(xy_data[[id_]]) & !is.null(rownames(xy_data))){
      rn = rownames(xy_data)
      xy_data = as.data.table(xy_data)
      xy_data[[id_]] = rn
    }else{
      xy_data = as.data.table(xy_data)
    }
  }
  if(is.data.frame(xy_data))
    if(is.null(xy_data[[set_]])){
      stop("set_ : '", set_, "' must be valid column in xy_data, not found!")
    }
  if(is.null(xy_data[[x_]])){
    stop("x_ : '", x_, "' must be valid column in xy_data, not found!")
  }
  if(is.null(xy_data[[y_]])){
    stop("y_ : '", y_, "' must be valid column in xy_data, not found!")
  }
  if(is.null(xy_data[[id_]])){
    stop("id_ : '", id_, "' must be valid column in xy_data, not found!")
  }
  if(is.null(xy_data[[set_density_]])){
    stop("set_density_ : '", set_density_, "' must be valid column in xy_data, not found!")
  }
  #labels
  if(is.na(main_title) || is.null(main_title)){
    main_title = ""
  }
  stopifnot(is.character(main_title))
  
  
  #set_ check
  if(!is.factor(xy_data[[set_]])){
    xy_data[[set_]] = factor(xy_data[[set_]])
  }
  sets.names = levels(xy_data[[set_]])
  sets.len = length(levels(xy_data[[set_]]))
  
  if(is.null(sets.colors)){
    sets.colors = RColorBrewer::brewer.pal(sets.len, "Dark2")
  }
  stopifnot(length(sets.colors) == sets.len)
  if(is.null(names(sets.colors))){
    names(sets.colors) = sets.names
  }
  if(bg.string %in% names(sets.colors) & is.character(bg.color)){
    sets.colors[bg.string] = bg.color
  }
  if(length(sets.sizes == 0)){
    sets.sizes = rep(sets.sizes, sets.len)
  }
  if(is.null(names(sets.sizes))){
    names(sets.sizes) = names(sets.colors)
  }
  stopifnot(names(sets.sizes) == sets.names)
  if(bg.string %in% names(sets.sizes) & is.numeric(bg.size)){
    sets.sizes[bg.string] = bg.size
  }
  
  # set_density_ check
  if(!is.factor(xy_data[[set_density_]])){
    xy_data[[set_density_]] = factor(xy_data[[set_density_]])
  }
  sets.density.names = levels(xy_data[[set_density_]])
  sets.density.len = length(levels(xy_data[[set_density_]]))
  
  if(is.null(sets.density.colors)){
    sets.density.colors = RColorBrewer::brewer.pal(sets.density.len, "Dark2")
  }
  stopifnot(length(sets.density.colors) == sets.density.len)
  if(is.null(names(sets.density.colors))){
    names(sets.density.colors) = sets.density.names
  }
  if(bg.density.string %in% names(sets.density.colors) & is.character(bg.density.color)){
    sets.density.colors[bg.density.string] = bg.density.color
  }
  if(length(sets.density.sizes == 0)){
    sets.density.sizes = rep(sets.density.sizes, sets.density.len)
  }
  if(is.null(names(sets.density.sizes))){
    names(sets.density.sizes) = names(sets.density.colors)
  }
  stopifnot(names(sets.density.sizes) == sets.density.names)
  if(bg.density.string %in% names(sets.density.sizes) & is.numeric(bg.density.size)){
    sets.sizes[bg.density.string] = bg.density.size
  }

  
  stopifnot(is.numeric(xy_data[[x_]]))
  stopifnot(is.numeric(xy_data[[y_]]))
  
  if(is.null(xlim_)){
    xlim = range(xy_data[[x_]])
  }else{
    xlim = xlim_
  }
  if(is.null(ylim_)){
    ylim = range(xy_data[[y_]])
  }else{
    ylim = ylim_
  }
  
  xy_data = xy_data[order(get(set_), decreasing = TRUE)]
  gene_o = xy_data[set_ != bg.string,][order(get(set_), decreasing = TRUE)][order(abs(get(x_) - get(y_)), decreasing = TRUE)][[id_]]
  if(n_auto_label > length(gene_o)) n_auto_label = length(gene_o)
  to_label = c(manual_label, gene_o[seq_len(n_auto_label)])
  
  p_scatter = ggplot(xy_data, aes_string(x = x_, y = y_, 
                                         color = set_, size = set_,
                                         label = id_)) + 
    geom_point() +
    scale_color_manual(values = sets.colors, drop = FALSE) +
    scale_size_manual(values = sets.sizes, drop = FALSE) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = labs_x, y = labs_y, color = labs_sets, size = labs_sets) +
    guides() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(t = -4, r = -4, b = 2, l = 2, unit = "pt"))
  p_x_density = ggplot(mapping = aes_string(x = x_, color = set_density_)) +
    geom_line(data = xy_data, stat = "density", color = bg.density.color) +
    geom_line(data = xy_data[get(set_density_) != bg.density.string], stat = "density") +
    scale_color_manual(values = sets.density.colors, drop = FALSE) +
    coord_cartesian(xlim = xlim) +
    labs(x = "", color = labs_sets_density) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  p_y_density = ggplot(mapping = aes_string(x = y_, color = set_density_)) +
    geom_line(data = xy_data, stat = "density", color = bg.density.color) +
    geom_line(data = xy_data[get(set_density_) != bg.density.string], stat = "density") +
    scale_color_manual(values = sets.density.colors, drop = FALSE) +
    coord_flip(xlim = ylim) +
    labs(x = "", color = labs_sets_density)+
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  
  #add reference lines
  if(is.numeric(ref_line.x)){
    if(length(ref_line.x.color) == 1){
      ref_line.x.color = rep(ref_line.x.color, length(ref_line.x))
    }
    stopifnot(length(ref_line.x.color) == length(ref_line.x))
    for(i in seq_along(ref_line.x)){
      p_scatter = p_scatter + annotate("line", x = rep(ref_line.x[i], 2), y = ylim, color = ref_line.x.color[i], size = .5, linetype = "dashed")
      p_x_density = p_x_density + annotate("line", x = rep(ref_line.x[i], 2), y = get_gg_yrange(p_x_density), color = ref_line.x.color[i], size = .5, linetype = "dashed")
    }
    
  }
  if(is.numeric(ref_line.y)){
    if(length(ref_line.y.color) == 1){
      ref_line.y.color = rep(ref_line.y.color, length(ref_line.y))
    }
    stopifnot(length(ref_line.y.color) == length(ref_line.y))
    for(i in seq_along(ref_line.y)){
      p_scatter = p_scatter + annotate("line", x = xlim, y = rep(ref_line.y[i], 2), color = ref_line.y.color[i], size = .5, linetype = "dashed")
      p_y_density = p_y_density + annotate("line", x = rep(ref_line.y[i], 2), y = get_gg_yrange(p_x_density), color = ref_line.x.color[i], size = .5, linetype = "dashed")    
    }
    
  }
  if(is.numeric(ref_line.slope)){
    if(max(xlim) > max(ylim)){
      right_pt = max(ylim)
    }else{
      right_pt = max(xlim)
    }
    if(min(xlim) < min(ylim)){
      left_pt = min(ylim)
    }else{
      left_pt = min(xlim)
    }
    p_scatter = p_scatter + annotate("line", x = c(left_pt, right_pt), y = c(left_pt, right_pt), color = ref_line.slope.color, size = .5, linetype = "dashed")
  }
  
  #add labels
  if(length(to_label) > 0){
    p_scatter = p_scatter + label_GEOM(data = xy_data[get(id_) %in% to_label], show.legend = FALSE, size = label_size)        
  }
  
  
  
  components = list(scatter = p_scatter, x_density = p_x_density, y_density = p_y_density)
  
  pg = plot_scatter_side_density_assemble(components, 
                                          main_title = main_title,
                                          main_title.x = main_title.x,
                                          main_title.y = main_title.y,
                                          main_title.hjust = main_title.hjust,
                                          main_title.vjust = main_title.vjust)
  
  if(!suppress_plot)
    plot(pg)
  invisible(list(assembled = pg, components = components))
}



# custom plot_scatter_side_density_xy ------------------------------------------
plot_scatter_side_density_xy_rel = function( xy_data,
                                         x_ = "x",
                                         y_ = "y",
                                         id_ = "id",
                                         set_ = "set",
                                         set_density_ = set_,
                                         #labels
                                         labs_x = x_,
                                         labs_y = y_,
                                         labs_sets = set_,
                                         labs_sets_density = set_density_,
                                         main_title = NULL,
                                         main_title.x = .02,
                                         main_title.y = .5,
                                         main_title.hjust = 0,
                                         main_title.vjust = .5,
                                         #point color and sizing
                                         sets.colors = NULL,
                                         sets.density.colors = NULL,
                                         bg.string = "background",
                                         bg.density.string = "background",
                                         bg.color = "gray70",
                                         bg.density.color = "gray70",
                                         sets.sizes = 1,
                                         sets.density.sizes = 1,
                                         bg.size = .5,
                                         bg.density.size = 0.5,
                                         #limits
                                         xlim_ = NULL,
                                         ylim_ = NULL,
                                         #point labelling
                                         n_auto_label = 8,
                                         manual_label = NULL,
                                         label_size = 6,
                                         label_color = 'black',
                                         label_GEOM = shadowtext::geom_shadowtext,
                                         #reference lines
                                         ref_line.x = 0,
                                         ref_line.x.color = "gray50",
                                         ref_line.y = 0,
                                         ref_line.y.color = "gray50",
                                         ref_line.slope = 1,
                                         ref_line.slope.color = "black",
                                         suppress_plot = FALSE){
  if(is.matrix(xy_data)){
    rn = rownames(xy_data)
    xy_data = as.data.table(xy_data)
    xy_data[[id_]] = rn
  }
  if(is.data.frame(xy_data) & !is.data.table(xy_data)){
    if(is.null(xy_data[[id_]]) & !is.null(rownames(xy_data))){
      rn = rownames(xy_data)
      xy_data = as.data.table(xy_data)
      xy_data[[id_]] = rn
    }else{
      xy_data = as.data.table(xy_data)
    }
  }
  if(is.data.frame(xy_data))
    if(is.null(xy_data[[set_]])){
      stop("set_ : '", set_, "' must be valid column in xy_data, not found!")
    }
  if(is.null(xy_data[[x_]])){
    stop("x_ : '", x_, "' must be valid column in xy_data, not found!")
  }
  if(is.null(xy_data[[y_]])){
    stop("y_ : '", y_, "' must be valid column in xy_data, not found!")
  }
  if(is.null(xy_data[[id_]])){
    stop("id_ : '", id_, "' must be valid column in xy_data, not found!")
  }
  if(is.null(xy_data[[set_density_]])){
    stop("set_density_ : '", set_density_, "' must be valid column in xy_data, not found!")
  }
  #labels
  if(is.na(main_title) || is.null(main_title)){
    main_title = ""
  }
  stopifnot(is.character(main_title))
  
  
  #set_ check
  if(!is.factor(xy_data[[set_]])){
    xy_data[[set_]] = factor(xy_data[[set_]])
  }
  sets.names = levels(xy_data[[set_]])
  sets.len = length(levels(xy_data[[set_]]))
  
  if(is.null(sets.colors)){
    sets.colors = RColorBrewer::brewer.pal(sets.len, "Dark2")
  }
  stopifnot(length(sets.colors) == sets.len)
  if(is.null(names(sets.colors))){
    names(sets.colors) = sets.names
  }
  if(bg.string %in% names(sets.colors) & is.character(bg.color)){
    sets.colors[bg.string] = bg.color
  }
  if(length(sets.sizes == 0)){
    sets.sizes = rep(sets.sizes, sets.len)
  }
  if(is.null(names(sets.sizes))){
    names(sets.sizes) = names(sets.colors)
  }
  stopifnot(names(sets.sizes) == sets.names)
  if(bg.string %in% names(sets.sizes) & is.numeric(bg.size)){
    sets.sizes[bg.string] = bg.size
  }
  
  # set_density_ check
  if(!is.factor(xy_data[[set_density_]])){
    xy_data[[set_density_]] = factor(xy_data[[set_density_]])
  }
  sets.density.names = levels(xy_data[[set_density_]])
  sets.density.len = length(levels(xy_data[[set_density_]]))
  
  if(is.null(sets.density.colors)){
    sets.density.colors = RColorBrewer::brewer.pal(sets.density.len, "Dark2")
  }
  stopifnot(length(sets.density.colors) == sets.density.len)
  if(is.null(names(sets.density.colors))){
    names(sets.density.colors) = sets.density.names
  }
  if(bg.density.string %in% names(sets.density.colors) & is.character(bg.density.color)){
    sets.density.colors[bg.density.string] = bg.density.color
  }
  if(length(sets.density.sizes == 0)){
    sets.density.sizes = rep(sets.density.sizes, sets.density.len)
  }
  if(is.null(names(sets.density.sizes))){
    names(sets.density.sizes) = names(sets.density.colors)
  }
  stopifnot(names(sets.density.sizes) == sets.density.names)
  if(bg.density.string %in% names(sets.density.sizes) & is.numeric(bg.density.size)){
    sets.sizes[bg.density.string] = bg.density.size
  }
  
  
  stopifnot(is.numeric(xy_data[[x_]]))
  stopifnot(is.numeric(xy_data[[y_]]))
  
  if(is.null(xlim_)){
    xlim = range(xy_data[[x_]])
  }else{
    xlim = xlim_
  }
  if(is.null(ylim_)){
    ylim = range(xy_data[[y_]])
  }else{
    ylim = ylim_
  }
  
  xy_data = xy_data[order(get(set_), decreasing = TRUE)]
  gene_o = xy_data[set_ != bg.string,][order(get(set_), decreasing = TRUE)][order(abs(get(x_) - get(y_)), decreasing = TRUE)][[id_]]
  if(n_auto_label > length(gene_o)) n_auto_label = length(gene_o)
  to_label = c(manual_label, gene_o[seq_len(n_auto_label)])
  
  p_scatter = ggplot(xy_data, aes_string(x = x_, y = y_, 
                                         color = set_, size = set_,
                                         label = id_)) + 
    geom_point() +
    scale_color_manual(values = sets.colors, drop = FALSE) +
    scale_size_manual(values = sets.sizes, drop = FALSE) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = labs_x, y = labs_y, color = labs_sets, size = labs_sets) +
    guides() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(t = -4, r = -4, b = 2, l = 2, unit = "pt"))
  p_x_density = ggplot(mapping = aes_string(x = x_, color = set_density_)) +
    geom_line(data = xy_data %>% dplyr::filter(gwas_g == "Other"), 
              stat = "density", color = bg.density.color) +
    geom_line(data = xy_data %>% dplyr::filter(gwas_g != "CTmax"), stat = "density") +
    scale_color_manual(values = sets.density.colors, drop = FALSE) +
    coord_cartesian(xlim = xlim) +
    labs(x = "", color = labs_sets_density) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  p_y_density = ggplot(mapping = aes_string(x = y_, color = set_density_)) +
    geom_line(data = xy_data %>% dplyr::filter(gwas_g == "Other"), 
              stat = "density", color = bg.density.color) +
    geom_line(data = xy_data %>% dplyr::filter(gwas_g != "CTmin"), stat = "density") +
    scale_color_manual(values = sets.density.colors, drop = FALSE) +
    coord_flip(xlim = ylim) +
    labs(x = "", color = labs_sets_density)+
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  
  #add reference lines
  if(is.numeric(ref_line.x)){
    if(length(ref_line.x.color) == 1){
      ref_line.x.color = rep(ref_line.x.color, length(ref_line.x))
    }
    stopifnot(length(ref_line.x.color) == length(ref_line.x))
    for(i in seq_along(ref_line.x)){
      p_scatter = p_scatter + annotate("line", x = rep(ref_line.x[i], 2), y = ylim, color = ref_line.x.color[i], size = .5, linetype = "dashed")
      p_x_density = p_x_density + annotate("line", x = rep(ref_line.x[i], 2), y = get_gg_yrange(p_x_density), color = ref_line.x.color[i], size = .5, linetype = "dashed")
    }
    
  }
  if(is.numeric(ref_line.y)){
    if(length(ref_line.y.color) == 1){
      ref_line.y.color = rep(ref_line.y.color, length(ref_line.y))
    }
    stopifnot(length(ref_line.y.color) == length(ref_line.y))
    for(i in seq_along(ref_line.y)){
      p_scatter = p_scatter + annotate("line", x = xlim, y = rep(ref_line.y[i], 2), color = ref_line.y.color[i], size = .5, linetype = "dashed")
      p_y_density = p_y_density + annotate("line", x = rep(ref_line.y[i], 2), y = get_gg_yrange(p_x_density), color = ref_line.x.color[i], size = .5, linetype = "dashed")    
    }
    
  }
  if(is.numeric(ref_line.slope)){
    if(max(xlim) > max(ylim)){
      right_pt = max(ylim)
    }else{
      right_pt = max(xlim)
    }
    if(min(xlim) < min(ylim)){
      left_pt = min(ylim)
    }else{
      left_pt = min(xlim)
    }
    p_scatter = p_scatter + annotate("line", x = c(left_pt, right_pt), y = c(left_pt, right_pt), color = ref_line.slope.color, size = .5, linetype = "dashed")
  }
  
  #add labels
  if(length(to_label) > 0){
    p_scatter = p_scatter + label_GEOM(data = xy_data[get(id_) %in% to_label], show.legend = FALSE, size = label_size)        
  }
  
  
  
  components = list(scatter = p_scatter, x_density = p_x_density, y_density = p_y_density)
  
  pg = plot_scatter_side_density_assemble(components, 
                                          main_title = main_title,
                                          main_title.x = main_title.x,
                                          main_title.y = main_title.y,
                                          main_title.hjust = main_title.hjust,
                                          main_title.vjust = main_title.vjust)
  
  if(!suppress_plot)
    plot(pg)
  invisible(list(assembled = pg, components = components))
}


# custom plot_scatter_side_density_assemble ------------------------------------
plot_scatter_side_density_assemble = function(components, 
                                              main_title = "", 
                                              main_title.x = .02, 
                                              main_title.y = .5, 
                                              main_title.hjust = 0, 
                                              main_title.vjust = .5){
  
  p_scatter = components$scatter
  p_x_density = components$x_density
  p_y_density = components$y_density
  
  
  p_legend = cowplot::get_legend(p_scatter)
  #d_legend = cowplot::get_legend(p_y_density)
  
  # legends = cowplot::plot_grid(plotlist = c(list(p_legend), list(d_legend)),
  #                              sscale = 0.5)
  
  grobs_y = sync_height(list(p_scatter + 
                               guides(color = "none", size = "none"), 
                             p_y_density + 
                               guides(color = "none")))
  
  grobs_x = sync_width(list(grobs_y[[1]], 
                            p_x_density + 
                              guides(color = "none")))
  
  
  
  pg = cowplot::plot_grid(plotlist = c(grobs_x[2], list(p_legend), grobs_y), 
                          rel_widths = c(2, 1), rel_heights = c(1,2))
  if(main_title != ""){
    pg = cowplot::plot_grid(
      cowplot::ggdraw() + 
        cowplot::draw_text(main_title, 
                           x = main_title.x, 
                           y = main_title.y, 
                           hjust = main_title.hjust, 
                           vjust = main_title.vjust),
      pg,
      rel_heights = c(1, 15),
      ncol = 1
    )
  }
  pg
}

# ggmaplot_custom --------------------------------------------------------------

# ggpubr functions that are needed
.levels <- function(x){
  if(!is.factor(x)) x <- as.factor(x)
  levels(x)
}

.parse_font <- function(font){
  if(is.null(font)) res <- NULL
  else if(inherits(font, "list")) res <- font
  else{
    # matching size and face
    size <- grep("^[0-9]+$", font, perl = TRUE)
    face <- grep("plain|bold|italic|bold.italic", font, perl = TRUE)
    if(length(size) == 0) size <- NULL else size <- as.numeric(font[size])
    if(length(face) == 0) face <- NULL else face <- font[face]
    color <- setdiff(font, c(size, face))
    if(length(color) == 0) color <- NULL
    res <- list(size=size, face = face, color = color)
  }
  res
}

ggmaplot_gwas <- function (data, fdr = 0.05, fc = 1.5, genenames = NULL,
                      detection_call = NULL, size = NULL,
                      font.label = c(12, "plain", "black"),
                      font.label.sig = c(12, "plain", "red"),
                      label.rectangle = FALSE,
                      palette = c("#B31B21", "#1465AC", "gray1", "darkgray"),
                      top = 15, select.top.method = c("padj", "fc", "gwas"),
                      main = NULL, xlab = "Log2 mean expression",  ylab = "Log2 fold change",
                      ggtheme = theme_classic(),...)
{
  
  if(!base::inherits(data, c("matrix", "data.frame", "DataFrame", "DE_Results", "DESeqResults")))
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  if(!is.null(detection_call)){
    if(nrow(data)!=length(detection_call))
      stop("detection_call must be a numeric vector of length = nrow(data)")
  }
  else if("detection_call" %in% colnames(data)){
    detection_call <- as.vector(data$detection_call)
  }
  else detection_call = rep(1, nrow(data))
  
  # Legend position
  if(is.null(list(...)$legend)) legend <- c(0.12, 0.9)
  
  # Check data format
  ss <- base::setdiff(c("baseMean", "log2FoldChange", "padj"), colnames(data))
  if(length(ss)>0) stop("The colnames of data must contain: ",
                        paste(ss, collapse = ", "))
  
  if(is.null(genenames)) genenames <- rownames(data)
  else if(length(genenames)!=nrow(data))
    stop("genenames should be of length nrow(data).")
  
  # Create the levels for coloring
  sig <- rep(3, nrow(data)) # all rows are defaulted to group 4
  
  sig[which(data$padj <= fdr & 
              data$log2FoldChange < 0 & 
              abs(data$log2FoldChange) >= log2(fc) & 
              detection_call == 1)] <- 2
  down_n <- sum(sig == 2)
  
  sig[which(data$padj <= fdr & 
              data$log2FoldChange > 0 & 
              abs(data$log2FoldChange) >= log2(fc) & 
              detection_call == 1)] <- 1
  up_n <- sum(sig == 1)
  
  sig[which(!is.na(data$ID))] <- 4
  gwas_n <- sum(sig == 4)
  
  # making a simpler data.frame with new names for each variable
  data <- data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange,
                     padj = data$padj, sig = sig)
  
  # Change level labels
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- .levels(data$sig) %>% as.numeric()
  #palette <- palette[.lev]
  new.levels <- c(
    paste0("Up: ", up_n),
    paste0("Down: ", down_n),
    "NS",
    paste0("GWAS: ", gwas_n)
  ) %>% .[.lev]
  
  data$sig <- factor(data$sig, labels = new.levels)
  
  
  # Ordering for selecting top gene
  select.top.method <- match.arg(select.top.method)
  if(select.top.method == "padj") data <- data[order(data$padj), ]
  else if(select.top.method == "fc") data <- data[order(abs(data$lfc), decreasing = TRUE), ]
  else if(select.top.method == "gwas") data <- data[order(data$sig, decreasing = TRUE), ]
  # select data for top genes
  labs_data <- stats::na.omit(data)
  #labs_data <- subset(labs_data, padj <= fdr & name!="" & abs(lfc) >= log2(fc))
  labs_data <- labs_data %>%
    dplyr::filter(sig == paste0("GWAS: ", gwas_n)) %>%
    dplyr::top_n(top)
  
  labs_data_sig <- labs_data %>%
    dplyr::filter(padj <= fdr)
  
  labs_data_ns <- labs_data %>%
    dplyr::filter(padj > fdr)
  
  
  font.label <- .parse_font(font.label)
  font.label$size <- ifelse(is.null(font.label$size), 12, font.label$size)
  font.label$color <- ifelse(is.null(font.label$color), "black", font.label$color)
  font.label$face <- ifelse(is.null(font.label$face), "plain", font.label$face)
  
  font.label.sig <- .parse_font(font.label.sig)
  font.label.sig$size <- ifelse(is.null(font.label.sig$size), 12, font.label.sig$size)
  font.label.sig$color <- ifelse(is.null(font.label.sig$color), "black", font.label.sig$color)
  font.label.sig$face <- ifelse(is.null(font.label.sig$face), "plain", font.label.sig$face)
  
  
  data <- data %>%
    dplyr::arrange(sig)
  
  # Plot
  set.seed(2)
  mean <- lfc <- sig <- name <- padj <-  NULL
  p <- ggplot(data = NULL, aes(x = log2(mean + 1), y = lfc)) +
    geom_point(data = data %>% 
                 dplyr::filter(sig != paste0("GWAS: ", gwas_n)), 
               aes(x = log2(mean + 1), y = lfc, color = sig), 
                size = size) + 
    geom_point(data = data %>% 
                 dplyr::filter(sig == paste0("GWAS: ", gwas_n)),
               aes(x = log2(mean + 1), y = lfc, color = sig),
               size = 2)
    
  if(label.rectangle){
    p <- p + 
      
      ggrepel::geom_label_repel(data = labs_data_sig, mapping = aes(label = name),
                                       box.padding = unit(0.35, "lines"),
                                       point.padding = unit(0.3, "lines"),
                                       force = 10, fontface = font.label.sig$face,
                                       size = font.label.sig$size/3, color = font.label.sig$color) +
                                #, nudge_x = -3, nudge_y = -2) +
      
      ggrepel::geom_label_repel(data = labs_data_ns, mapping = aes(label = name),
                                box.padding = unit(0.35, "lines"),
                                point.padding = unit(0.3, "lines"),
                                force = 10, fontface = font.label$face,
                                size = font.label$size/3, color = font.label$color) 
                                #, nudge_x = 3)
  }
  else{
    p <- p + 
      
      ggrepel::geom_text_repel(data = labs_data_sig, mapping = aes(label = name),
                                      box.padding = unit(0.35, "lines"),
                                      point.padding = unit(0.3, "lines"),
                                      force = 10, fontface = font.label$face,
                                      size = font.label$size/3, color = "red") +
      
      ggrepel::geom_text_repel(data = labs_data_ns, mapping = aes(label = name),
                               box.padding = unit(0.35, "lines"),
                               point.padding = unit(0.3, "lines"),
                               force = 10, fontface = font.label.sig$face,
                               size = font.label.sig$size/3, color = font.label.sig$color)
  }
  
  p <- p + scale_x_continuous(breaks = seq(0, max(log2(data$mean + 1)), 2))+
    labs(x = xlab, y = ylab, title = main, color = "")+ # to remove legend title use color = ""
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2),
               color = c("black", "black", "black")) 
  
  p <- ggpar(p, palette = palette, ggtheme = ggtheme, ...) 
  p + scale_color_manual(breaks = new.levels, values = palette)
}




#################TEST



ggmaplot_test <- function (data, fdr = 0.05, fc = 1.5, genenames = NULL,
                           detection_call = NULL, size = NULL,
                           font.label = c(12, "plain", "black"),
                           font.label.sig = c(12, "plain", "red"),
                           label.rectangle = FALSE,
                           palette = c("#B31B21", "#1465AC", "gray1", "darkgray"),
                           color_sig = FALSE,
                           top = 15, select.top.method = c("padj", "fc", "gwas"),
                           main = NULL, xlab = "Log2 mean expression",  ylab = "Log2 fold change",
                           ggtheme = theme_classic(),...)
{
  
  if(!base::inherits(data, c("matrix", "data.frame", "DataFrame", "DE_Results", "DESeqResults")))
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  if(!is.null(detection_call)){
    if(nrow(data)!=length(detection_call))
      stop("detection_call must be a numeric vector of length = nrow(data)")
  }
  else if("detection_call" %in% colnames(data)){
    detection_call <- as.vector(data$detection_call)
  }
  else detection_call = rep(1, nrow(data))
  
  # Legend position
  if(is.null(list(...)$legend)) legend <- c(0.12, 0.9)
  
  # Check data format
  ss <- base::setdiff(c("baseMean", "log2FoldChange", "padj"), colnames(data))
  if(length(ss)>0) stop("The colnames of data must contain: ",
                        paste(ss, collapse = ", "))
  
  if(is.null(genenames)) genenames <- rownames(data)
  else if(length(genenames)!=nrow(data))
    stop("genenames should be of length nrow(data).")
  
  # Create the levels for coloring
  sig <- rep(3, nrow(data)) # all rows are defaulted to group 4
  
  sig[which(data$padj <= fdr & 
              data$log2FoldChange < 0 & 
              abs(data$log2FoldChange) >= log2(fc) & 
              detection_call == 1)] <- 2
  down_n <- sum(sig == 2)
  
  sig[which(data$padj <= fdr & 
              data$log2FoldChange > 0 & 
              abs(data$log2FoldChange) >= log2(fc) & 
              detection_call == 1)] <- 1
  up_n <- sum(sig == 1)
  
  sig[which(!is.na(data$ID))] <- 4
  gwas_n <- sum(sig == 4)
  
  # making a simpler data.frame with new names for each variable
  data <- data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange,
                     padj = data$padj, sig = sig)
  
  # Change level labels
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- .levels(data$sig) %>% as.numeric()
  palette <- palette[.lev]
  new.levels <- c(
    paste0("Up: ", up_n),
    paste0("Down: ", down_n),
    "NS",
    paste0("GWAS: ", gwas_n)
  ) %>% .[.lev]
  
  data$sig <- factor(data$sig, labels = new.levels)
  
  
  # Ordering for selecting top gene
  select.top.method <- match.arg(select.top.method)
  if(select.top.method == "padj") data <- data[order(data$padj), ]
  else if(select.top.method == "fc") data <- data[order(abs(data$lfc), decreasing = TRUE), ]
  else if(select.top.method == "gwas") data <- data[order(data$sig, decreasing = TRUE), ]
  # select data for top genes
  labs_data <- stats::na.omit(data)
  #labs_data <- subset(labs_data, padj <= fdr & name!="" & abs(lfc) >= log2(fc))
  labs_data <- labs_data %>%
    dplyr::filter(sig == paste0("GWAS: ", gwas_n)) %>%
    dplyr::top_n(top)
  
  font.label <- .parse_font(font.label)
  font.label$size <- ifelse(is.null(font.label$size), 12, font.label$size)
  font.label$color <- ifelse(is.null(font.label$color), "black", font.label$color)
  font.label$face <- ifelse(is.null(font.label$face), "plain", font.label$face)
  
  font.label.sig <- .parse_font(font.label.sig)
  font.label.sig$size <- ifelse(is.null(font.label.sig$size), 12, font.label.sig$size)
  font.label.sig$color <- ifelse(is.null(font.label.sig$color), "black", font.label.sig$color)
  font.label.sig$face <- ifelse(is.null(font.label.sig$face), "plain", font.label.sig$face)
  
  
  data <- data %>%
    dplyr::arrange(sig)
  
  # Plot
  set.seed(2)
  mean <- lfc <- sig <- name <- padj <-  NULL
  p <- ggplot(data = NULL, aes(x = log2(mean + 1), y = lfc)) +
    geom_point(data = data %>% 
                 dplyr::filter(sig != paste0("GWAS: ", gwas_n)), 
               aes(x = log2(mean + 1), y = lfc, color = sig), 
               size = size) + 
    geom_point(data = data %>% 
                 dplyr::filter(sig == paste0("GWAS: ", gwas_n)),
               aes(x = log2(mean + 1), y = lfc, color = sig),
               size = 2)
  
  
  p <- p + scale_x_continuous(breaks = seq(0, max(log2(data$mean + 1)), 2))+
    labs(x = xlab, y = ylab, title = main, color = "")+ # to remove legend title use color = ""
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2),
               color = c("black", "black", "black")) 
  
  
  if(label.rectangle){
    p <- p + 
      
      ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = name, color = labs_data$padj < fdr),
                                box.padding = unit(0.35, "lines"),
                                point.padding = unit(0.3, "lines"),
                                force = 10, fontface = font.label.sig$face,
                                size = font.label.sig$size/3) +
      guides(colour = guide_legend(override.aes = list(shape = 20)))
      }
  
  
  else{
    p <- p + 
      
      ggrepel::geom_text_repel(data = labs_data_sig, mapping = aes(label = name, color = labs_data$padj < fdr),
                               box.padding = unit(0.35, "lines"),
                               point.padding = unit(0.3, "lines"),
                               force = 10, fontface = font.label$face,
                               size = font.label$size/3, color = "red")
  }
  
  p <- ggpar(p, ggtheme = ggtheme, ...) 
  
  if (color_sig) {
    p <- p + scale_color_manual(breaks = new.levels, values = c("#1465AC", "black", "black", "#A6A6A680", "firebrick", "#B31B21"))
  } 
  else {
    p <- p + scale_color_manual(breaks = new.levels, values = c("#1465AC", "black", "black", "#A6A6A680", "black", "#B31B21"))
  }
  #p + scale_color_manual(breaks = new.levels, values = c("#1465AC", "black", "black", "#A6A6A680", "firebrick", "#B31B21"))
  #p + scale_color_manual(breaks = new.levels, values = palette)
  
  p
}


# function to get the gene from the dgrp .top.annot file -----------------------

gwas_dgrp_annot <- function(dat){
  # extract only the bit at the beginning
  x <- stringr::str_extract(dat$GeneAnnotation, 
                            "SiteClass\\[[\\w\\W|]+\\],")
  
  # remove the SiteClass and brackets etc.
  dat$annot <- stringr::str_remove(stringr::str_remove(x, "SiteClass\\["), 
                                   "\\],")
  
  # split up different genes and save fbgn, gene, feature, and pos in new columns
  df <- dat %>%
    tidyr::separate(col = annot, into = c("gene1", "gene2"), sep = ";") %>%
    tidyr::pivot_longer(gene1:gene2, names_to = "lab", values_to = "gene",
                        values_drop_na = TRUE) %>%
    tidyr::separate(col = gene, 
                    into = c("fbgn", "gene", "feature", "pos"), sep = "\\|") 
  
  return(df)
}

# function to get the gene from the dgrp .top.annot file -----------------------

gwas_dgrp_reg_annot <- function(dat){
  
  # remove the parens
  dat$annot <- stringr::str_remove_all(stringr::str_remove_all(dat$RegulationAnnotation, "\\("),
                                       "\\)")
  # get the max number of regulatory annotations to make column names when split 
  num_regs <- max(stringr::str_count(dat$annot, pattern = ";") + 1)
  
  # split up different genes and save fbgn, gene, feature, and pos in new columns
  df <- dat %>%
    tidyr::separate(col = annot, into = paste("annot", 1:num_regs, sep = "_"), sep = ";") %>%
    tidyr::pivot_longer(dplyr::contains("annot_"), names_to = "lab_reg", values_to = "reg",
                        values_drop_na = TRUE) %>%
    tidyr::separate(col = reg, 
                    into = c("reg_type", "reg_source", "reg_flybaseID"), sep = "\\|") 
  
  return(df)
}


# function to group the genes by differential expression -----------------------

group_deg <- function(dat, pval_cut = 0.01, lfc_cut = 0){
  
  deg_grouped <- dat %>%
    mutate(cold_genes = case_when(padj.cold < pval_cut ~ "sig",
                                  padj.cold >= pval_cut ~ "ns"),
           hot_genes = case_when(padj.hot < pval_cut ~ "sig",
                                 padj.hot >= pval_cut ~ "ns"),
           cold_fc = case_when(log2FoldChange.cold < lfc_cut ~ "down",
                               log2FoldChange.cold > lfc_cut ~ "up"),
           hot_fc = case_when(log2FoldChange.hot < lfc_cut ~ "down",
                              log2FoldChange.hot > lfc_cut ~ "up")) %>%
    filter(!is.na(cold_genes) & !is.na(hot_genes) &
             !is.na(cold_fc) & !is.na(hot_fc)) %>%
    mutate(col = paste(cold_genes, hot_genes, cold_fc, hot_fc)) %>%
    mutate(group = case_when(col == "sig sig up up" |
                               col == "sig sig down down" ~ "Shared",
                             col == "sig ns up up" |
                               col == "ns sig up up" |
                               col == "sig ns down down" |
                               col == "ns sig down down" ~ "Sig1 Shared",
                             col == "ns ns down down" |
                               col == "ns ns down up" |
                               col == "ns ns up down" |
                               col == "ns ns up up" ~ "NS",
                             col == "sig sig down up" |
                               col == "sig sig up down" ~ "Unique",
                             col == "ns sig up down" |
                               col == "ns sig down up" |
                               col == "sig ns down up" |
                               col == "sig ns up down" ~ "Sig1 Unique")) %>%
    add_count(group) %>% 
    mutate(groupn = paste0(group, ': ', n)) %>%
    arrange(groupn)
  
  return(deg_grouped)
}


# funcion to produce the data frame with the ks pval pairwise comparisons ------ 
ks.pairwise <- function(.data, dist) {
  x <- .data %>%
    dplyr::select(dist, gwas_g) %>%
    dplyr::group_by(gwas_g) %>%
    tidyr::pivot_wider(names_from = gwas_g, 
                       values_from = dist)
  
  cols <- colnames(x)
  
  comparisons <- combn(1:length(cols), 2)
  
  df <- NULL
  
  for (i in 1:ncol(comparisons)) {
    # pval of each comparison
    pval <- ks.test(unlist(x[, comparisons[1, i]]), 
                    unlist(x[, comparisons[2, i]]))$p.value
    # character string of the relevant comparison
    Comp <- paste(cols[comparisons[1, i]], "vs.", cols[comparisons[2, i]]) 
    # bind together
    df <- as.data.frame(rbind(df, cbind(Comp, pval)))
  }
  
  # treat pvals as numeric
  df$pval <- as.numeric(as.character(df$pval))
  
  return(df)
}


# function to get gobp results into long format --------------------------------
go_gene_match <- function(dat, dat_match, dat_deg, pval_cut = 10^-4){
  
  # dat_match pval cut off
  dat_match <- dat_match %>%
    dplyr::filter(AvgMixedPval < pval_cut) %>%
    dplyr::arrange(AvgMixedPval) %>%
    dplyr::distinct(.keep_all = TRUE)
  
  # get the max number of regulatory annotations to make column names when split 
  num_regs <- max(stringr::str_count(dat$genes, pattern = ";") + 1)
  
  # split up different genes and save fbgn, gene, feature, and pos in new columns
  df <- dat %>%
    tidyr::separate(col = genes, into = paste("gene", 1:num_regs, sep = "_"), sep = ";") %>%
    tidyr::pivot_longer(dplyr::contains("gene_"), names_to = "lab", values_to = "gene",
                        values_drop_na = TRUE) %>%
    dplyr::filter(gene %in% dat_match$gene) %>%
    dplyr::left_join(dat_match, by = "gene") %>%
    dplyr::select(GO, gene, AvgMixedPval) %>%
    dplyr::left_join(dat_deg, by = "gene") %>%
    dplyr::select(GO, gene, AvgMixedPval, log2FoldChange, padj) %>%
    dplyr::mutate(gene = paste(gene, AvgMixedPval, log2FoldChange, padj, sep = "_")) %>%
    dplyr::select(GO, gene) %>%
    tidyr::pivot_wider(names_from = GO, 
                       values_from = gene,
                       values_fn = list(gene = list)) %>%
    tidyr::pivot_longer(everything(),
                        names_to = "GO",
                        values_to = "gene_list")
  
  # make the gene_list a list separated by ;
  df$genes <- vector(mode = "character", length = nrow(df))
  for (i in 1:nrow(df)){
    df$genes[i] <- paste(unlist(df$gene_list[i]), sep = "", collapse = ";")
  }
  
  df <- df %>%
    dplyr::select(GO, genes) 

  return(df) 
}

# match pvals, lfc, and go match -----------------------------------------------

match_p_lfc <- function(dat, dat_ct, dat_deg){
  # keep only the snp with the lowest pval for each gene
  dat_ct <- dat_ct %>%
    arrange(AvgMixedPval) %>%
    distinct(gene, .keep_all = TRUE)
  
  # split up different genes and match with corresponding pvals, lfc, etc.
  df <- dat %>%
    tidyr::separate(col = gene, into = paste("gene", 1:num_regs, sep = "_"), sep = ";") %>%
    tidyr::pivot_longer(dplyr::contains("gene_"), names_to = "lab", values_to = "gene",
                        values_drop_na = TRUE) %>%
    dplyr::left_join(dat_ct, by = "gene") %>%
    dplyr::select(GO, gene, AvgMixedPval) %>%
    dplyr::left_join(dat_deg, by = "gene") %>%
    dplyr::select(GO, gene, AvgMixedPval, log2FoldChange, padj) %>%
    dplyr::arrange(AvgMixedPval) %>%
    dplyr::mutate(gene = paste(gene, AvgMixedPval, log2FoldChange, padj, sep = "_")) %>%
    dplyr::select(GO, gene) %>%
    tidyr::pivot_wider(names_from = GO, 
                       values_from = gene,
                       values_fn = list(gene = list)) %>%
    tidyr::pivot_longer(everything(),
                        names_to = "GO",
                        values_to = "gene_list")
  
  # make the gene_list a list separated by ;
  df$genes <- vector(mode = "character", length = nrow(df))
  for (i in 1:nrow(df)){
    df$genes[i] <- paste(unlist(df$gene_list[i]), sep = "", collapse = ";")
  }
  
  df <- df %>%
    dplyr::select(GO, genes) %>%
    tidyr::separate(col = GO, into = c("GO", "match", "FDR"), sep = ";") %>%
    dplyr::arrange(FDR)
  
  return(df) 
}



# ------------------------------------------------------------------------------
# Function: read_clean_gtf
# Description: Read in and clean gtf file
# Inputs: gtf_file
# Outputs: gtf cleaned data frame

read_clean_gtf <- function(file) {
  # Read in the gtf
  x <- read_delim(file, 
                  delim = "\t", 
                  col_names = c("chr", "source", "feature", "start", "end", 
                                "score", "strand", "frame", "attributes"))
  
  # Clean up and separate the attributes file
  x$attributes <- str_replace_all(x$attributes, "\"", "")
  x <- x %>%
    separate(attributes, 
             sep = ";", 
             into = c("gene_id", "gene_symbol", "transcript_id", 
                      "transcript_symbol", "notes", "extra", "extraextra"))
  # remove the name and the awkward space before
  x$gene_id <- str_replace_all(x$gene_id, 
                               "gene_id ", "")
  x$gene_symbol <- str_replace_all(x$gene_symbol, 
                                   " gene_symbol ", "")
  x$transcript_id <- str_replace_all(x$transcript_id, 
                                     " transcript_id ", "")
  x$transcript_symbol <- str_replace_all(x$transcript_symbol, 
                                         " transcript_symbol ", "")
  
  return(x)
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: transcript_to_gene
# Description: convert a vector of transcripts to their gene symbols
# Inputs: character vector of transcripts and gtf data frame
# Outputs: gene

require(tidyverse)

transcript_to_gene <- function (trans_vec, gtf_dat){
  gtf_dat <- gtf_dat %>%
    distinct(transcript_id, .keep_all = TRUE) %>%
    filter(transcript_id != "")
  dat <- tibble::enframe(trans_vec, name = NULL, value = "transcript_id")
  dat$gene <- dat$transcript_id
  for (i in 1:nrow(dat)){
    print(i)
    temp <- which(dat$transcript_id[i] == gtf_dat$transcript_id, TRUE)
    if (length(temp) == 0) {
      gene_replace <- NA
    } else {
      gene_replace <- gtf_dat$gene_symbol[temp]
    }
    dat$gene <- gsub(dat$transcript_id[i],
                     gene_replace,
                     dat$gene)
  }
  dat_wo_na <- dat$gene[which(!is.na(dat$gene))]
  return(dat_wo_na)
}
# End function -----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Function: transcript_to_gene_df
# Description: convert a vector of transcripts to their gene symbols
# Inputs: character vector of transcripts and gtf data frame
# Outputs: gene

require(tidyverse)

transcript_to_gene_df <- function (dat, gtf_dat){
  gtf_dat <- gtf_dat %>%
    distinct(transcript_id, .keep_all = TRUE) %>%
    filter(transcript_id != "")

  dat$gene <- dat$transcript_id
  for (i in 1:nrow(dat)){
    print(i)
    temp <- which(dat$transcript_id[i] == gtf_dat$transcript_id, TRUE)
    if (length(temp) == 0) {
      gene_replace <- NA
    } else {
      gene_replace <- gtf_dat$gene_symbol[temp]
    }
    dat$gene <- gsub(dat$transcript_id[i],
                     gene_replace,
                     dat$gene)
  }
  return(dat)
}
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: gene_assoc_window
# Description: Annotate with all genes associated with the specific window
# Inputs: fst data frame and gtf data frame
# Outputs: fst data frame with new gene_assoc column

require(tidyverse)

gene_assoc_window <- function(dat, gtf_dat) {
  
  dat$gene_assoc <- vector(mode = "character", length = nrow(dat))
  
  gtf_dat <- gtf_dat %>%
    filter(feature == "gene")
  
  
  for (i in 1:nrow(dat)){
    print(i)
    
    # Find all genes that overlap the window in anyway
    temp <- c( 
      which(dat$CHR[i] == gtf_dat$chr & 
              gtf_dat$end >= dat$start[i] & 
              gtf_dat$end <= dat$end[i], TRUE),
      which(dat$CHR[i] == gtf_dat$chr & 
              gtf_dat$start >= dat$start[i] & 
              gtf_dat$start <= dat$end[i], TRUE),
      which(dat$CHR[i] == gtf_dat$chr & 
              gtf_dat$start <= dat$start[i] & 
              gtf_dat$end >= dat$end[i], TRUE)
    )
    
    # Paste all gene symbols together with ; between. NA if there is none
    if (is_empty(temp)){
      dat$gene_assoc[i] <- NA
    } else{
      dat$gene_assoc[i] <- paste(unique(gtf_dat$gene_symbol[temp]), 
                                 collapse = ";")
    }
    
  }
  return(dat)
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: gene_assoc_snp
# Description: Annotate with all genes associated with the snp
# Inputs: fst data frame and gtf data frame
# Outputs: fst data frame with new gene_assoc column

require(tidyverse)

gene_assoc_snp <- function(dat, gtf_dat) {
  
  # Initialize a column to store the gene associations
  dat$gene_assoc <- vector(mode = "character", length = nrow(dat))
  
  # Loop through each row of the dataframe individually
  pb <- progress::progress_bar$new(total = nrow(dat))
  
  for (i in 1:nrow(dat)){
    # Print the row just to watch the progress
    pb$tick()
    
    # Find all genes that overlap the window in anyway
    temp <- c( 
      which(dat$CHR[i] == gtf_dat$chr &  
              gtf_dat$start <= dat$BP[i] & 
              gtf_dat$end >= dat$BP[i], TRUE)
    )
    
    # Paste all gene symbols together with ; between. NA if there is none
    if (is_empty(temp)){
      dat$gene_assoc[i] <- NA
    } else{
      x <- gtf_dat[temp, ]
      
      # Split the df based on gene_symbol
      y <- x %>% 
        group_split(gene_symbol)
      
      # Initialize an empty object to store the string
      genes_feature_str <- NULL
      
      # Loop through each gene
      for (j in length(y)){
        # Get the gene and feature information. Features separated by a comma
        gene_str <- unique(y[[j]]$gene_symbol)
        feature_str <- paste(unique(y[[j]]$feature), 
                             collapse = ",") # maybe add something to get the heirachy of features later
        # Paste the gene symbol to the associated features separated by a semicolon
        gene_feature_str <- paste(gene_str, feature_str, sep = ";")
        
        # Paste the genes together separated by a 
        genes_feature_str <- c(genes_feature_str, gene_feature_str)
      }
      
      dat$gene_assoc[i] <- paste(genes_feature_str, 
                                 collapse = "|")
      
    }
    
  }
  return(dat)
} 
# End function -----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Function: get_feature_info
# Description: description
# Inputs: input_description
# Outputs: output_description

require(tidyverse)

get_feature_info <- function(chr, bp){
  
  # Find the rows that are within the start - end range of each feature
  temp <- c( 
    which(chr == gtf_dat$chr &  
            gtf_dat$start <= bp & 
            gtf_dat$end >= bp, TRUE)
  )
  
  # Paste all gene symbols together with ; between. NA if there is none
  if (is_empty(temp)){
    features <- NA
  } else{
    x <- gtf_dat[temp, ]
    
    # Split the df based on gene_symbol
    y <- x %>% 
      group_split(gene_symbol)
    
    # Initialize an empty object to store the string
    genes_feature_str <- NULL
    
    # Loop through each gene
    for (j in length(y)){
      # Get the gene and feature information. Features separated by a comma
      gene_str <- unique(y[[j]]$gene_symbol)
      feature_str <- paste(unique(y[[j]]$feature), 
                           collapse = ",") # maybe add something to get the heirachy of features later
      # Paste the gene symbol to the associated features separated by a semicolon
      gene_feature_str <- paste(gene_str, feature_str, sep = ";")
      
      # Paste the genes together separated by a 
      genes_feature_str <- c(genes_feature_str, gene_feature_str)
      
    }
    
    
    features <- paste(genes_feature_str, collapse = "|")
    
  }
  return(features)
}
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: gene_assoc_snp_pmap
# Description: Annotate with all genes associated with the snp
# Inputs: fst data frame and gtf data frame
# Outputs: fst data frame with new gene_assoc column

require(tidyverse)

gene_assoc_snp_pmap <- function(dat, gtf_dat) {
    
  dat$gene_assoc <- pmap_chr(list(dat$CHR, dat$BP), get_feature_info)

  return(dat)
} 
# End function -----------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: gtf_up_down_stream_annot
# Description: description
# Inputs: input_description
# Outputs: output_description

gtf_up_down_stream_annot <- function(gtf_df, bp_up = 1000, bp_down = bp_up) {
  
  gtf_up_pos <- gtf_df %>%
    filter(feature == "gene", strand == "+") %>%
    mutate(end = start - 1,
           start = start - bp_up - 1,
           feature = "upstream")
  
  gtf_down_pos <- gtf_df %>%
    filter(feature == "gene", strand == "+") %>%
    mutate(start = end + 1,
           end = end + bp_down + 1,
           feature = "downstream")
  
  gtf_up_neg <- gtf_df %>%
    filter(feature == "gene", strand == "-") %>%
    mutate(end = start - 1,
           start = start - bp_down - 1,
           feature = "downstream")
  
  gtf_down_neg <- gtf_df %>%
    filter(feature == "gene", strand == "-") %>%
    mutate(start = end + 1,
           end = end + bp_up + 1,
           feature = "upstream")
  
  gtf_df <- bind_rows(gtf_down_pos, 
                      gtf_up_pos,
                      gtf_down_neg, 
                      gtf_up_neg,
                      gtf_df) %>% 
    arrange(desc(chr), start)
  
  return(gtf_df)
} 

# End function -----------------------------------------------------------------



