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


# plot gwas density on scatter plot

plot_scatter_side_density.xy.TSO = function( xy_data,
                                         x_ = "x",
                                         y_ = "y",
                                         id_ = "id",
                                         set_ = "set",
                                         set_density = "set",
                                         #labels
                                         labs_x = x_,
                                         labs_y = y_,
                                         labs_sets = set_,
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
  if(is.null(xy_data[[set_density]])){
    stop("set_density : '", set_density, "' must be valid column in xy_data, not found!")
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
  
  # set_density check
  if(!is.factor(xy_data[[set_density]])){
    xy_data[[set_density]] = factor(xy_data[[set_density]])
  }
  sets.density.names = levels(xy_data[[set_density]])
  sets.density.len = length(levels(xy_data[[set_density]]))
  
  if(is.null(sets.density.colors)){
    sets.density.colors = RColorBrewer::brewer.pal(sets.density.len, "Dark2")
  }
  stopifnot(length(sets.density.colors) == sets.density.len)
  if(is.null(names(sets.density.colors))){
    names(sets.density.colors) = sets.density.names
  }
  if(bg.density.string %in% names(sets.density.colors) & is.character(bg.density.color)){
    sets.colors[bg.density.string] = bg.density.color
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
  gene_o = xy_data[set != bg.string,][order(get(set_), decreasing = TRUE)][order(abs(get(x_) - get(y_)), decreasing = TRUE)][[id_]]
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
    theme(plot.margin = margin(t = -4, r = -4, b = 0, l = 0, unit = "pt"))
  p_x_density = ggplot(mapping = aes_string(x = x_, color = set_density)) +
    geom_density(data = xy_data, color = bg.density.color) +
    geom_density(data = xy_data[get(set_density) != bg.density.string]) +
    geom_hline(yintercept = 0, colour = "white", size = 1) +
    scale_color_manual(values = sets.density.colors, drop = FALSE) +
    coord_cartesian(xlim = xlim) +
    labs(x = "")+
    theme_classic() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  p_y_density = ggplot(mapping = aes_string(x = y_, color = set_density)) +
    geom_density(data = xy_data, color = bg.density.color) +
    geom_density(data = xy_data[get(set_density) != bg.density.string]) +
    geom_hline(yintercept = 0, colour = "white", size = 1) +
    scale_color_manual(values = sets.density.colors, drop = FALSE) +
    coord_flip(xlim = ylim) +
    labs(x = "")+
    theme_classic() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  
  # # grob the panel
  # p_x_den <- ggplotGrob(p_x_d)
  # p_x_density <- gtable::gtable_filter(p_x_den, "panel")
  
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
  
  pg = plot_scatter_side_density.assemble(components, 
                                          main_title = main_title,
                                          main_title.x = main_title.x,
                                          main_title.y = main_title.y,
                                          main_title.hjust = main_title.hjust,
                                          main_title.vjust = main_title.vjust)
  
  if(!suppress_plot)
    plot(pg)
  invisible(list(assembled = pg, components = components))
}
