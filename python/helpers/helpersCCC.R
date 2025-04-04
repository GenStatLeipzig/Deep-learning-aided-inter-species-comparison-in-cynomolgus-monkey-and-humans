############### Module containing helper functions for CCC noteboks ###############
# Some functions were taken from other sources and modified as needed             #
# In these cases it is noted in the function where the originial was taken from   #
# Author: kmlr81, Kristina Müller                                                 #
###################################################################################


#' Generate a circus plot displaying the number of CCC events between different celltypes for a condition/timepoint in a dataset as detected by liana's rank_aggregate() function. 
#' Author: Holger Kristen, edited by: Kristina Müller
#' 
#' @param data A data frame with the columns 'source', 'target' and 'number' for, respectively, the source celltype, the target celltype and the number of detected CCC interactions based on ligand/target expression.
#' @param plot_title String, title for the generated circos plot. Default = "".
#' @param color_fill_arrows String vector containing two ggplot2 recognized color names. The colors indicating number of interactions for one connection between celltypes scaled from low to high using these two colors as min and max points. Default = c("yellow", "purple").
#' @param cluster_colors A string vector containing a number of ggplot2 recognized color names equal to the number of celltypes/clusters in the data. Assigns a given color to each celltype. Default = NULL.
#' @param legendsize A number. Sets the legend size of the plot accordingly. Default = 1.
#' @param min_range A number (float). Sets the min number of connections. Can be used to set a uniform colorscale for the number of CCC interactions between celltypes across multiple plots. Default = NULL. 
#' @param max_range A number (float). Sets the max number of connections. Can be used to set a uniform colorscale for the number of CCC interactions between celltypes across multiple plots. Default = NULL. 
#' @examples
#' vizualize_interactions(counts_df, plot_title="Circos plot for condition healthy", color_fill_arrows=c("blue", "green"), cluster_colors=c("orange", "red", "gray", "purple"), legendsize=2, min_reange=0, max_range=300)
visualize_interactions <- function(data,
                                   plot_title = "",
                                   color_fill_arrows = c("yellow", "purple"),
                                   cluster_colors = NULL, 
                                   legendsize= 1,
                                   min_range= NULL,
                                   max_range= NULL) {

  # modified from https://github.com/SCA-IRCM/SingleCellSignalR/blob/master/R/visualize.R with some help from https://claude.ai/ HK 27.5.24
  # Convert the input data frame to a data.table
  setDT(data)
  
  # Create a vector of unique cluster names
  clusters <- unique(c(data$source_cluster, data$target_cluster))
  
  # Assign default colors to each cluster if cluster_colors is not provided
  if (is.null(cluster_colors)) {
    cols <- rainbow(length(clusters))
    names(cols) <- clusters
  } else {
    if (length(cluster_colors) < length(clusters)) {
      stop("The number of cluster colors provided must be at least the number of clusters, i.e. ",length(clusters))
    }
    cols <- cluster_colors
    names(cols) <- clusters
  }
  
  # Prepare the data for the chord diagram
  dat <- data[, .(value = sum(number)), by = .(source_cluster, target_cluster)]
  
  # Create a color palette based on the range of interaction numbers or given min/max values
  if (is.null(min_range)) {
      min_range = min(dat$value)
  }
  if (is.null(max_range)) {
      max_range = max(dat$value)
  }
    
  color_palette <- colorRamp2(c(min_range, max_range), color_fill_arrows)
  
  # Create the chord diagram
  chordDiagramFromDataFrame(dat, annotationTrack = "grid", preAllocateTracks = 1,
                            directional = 1, direction.type = "arrows",
                            link.arr.length = 0.3, link.arr.width = 0.3,
                            link.arr.type = "triangle", link.arr.lty = par("lty"),
                            link.arr.lwd = 2.5, link.arr.col = alpha("black", 0.6),
                            grid.col = cols, col = color_palette(dat$value), big.gap = 10, small.gap = 2)
  
  # Add labels and axis to the chord diagram
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + 0.1, sector.name, facing = "bending.inside",
                niceFacing = FALSE, adj = c(0.5, 0), cex = 1.1)
    circos.axis(h = "top", labels = FALSE, major.tick = FALSE,
                major.tick.length = 1, sector.index = sector.name,
                track.index = 2)
  }, bg.border = NA)
  
  # Add a legend for the interaction numbers
  color_legend <- color_palette(seq(min_range, max_range, length.out = 5))
  legend_labels <- round(seq(min_range, max_range, length.out = 5))
  legend(x = "topleft", legend = legend_labels, fill = color_legend, title = "Number of \nInteractions", cex = legendsize )
  title(main=list(plot_title, cex = 1.3), line=-1.5)

}



#' Generate a vertical diverging barplot to show if the number of inferred CCC interactions between celltype pairs has increased or decreased (and by how much) between two conditions/timepoints in a dataset.
#' Author: Kristina Müller (kmlr81)
#'
#' @param data A dataframe with three columns: "source", "target" and "number". Containing respectively the name of the source cluster of the CCC interaction, the name of the target cluster and the difference in detected interactions between the two conditions/timepoints.
#' @param title A string to be set as the title of the plot. Default = "".
#' @param colors A vector containing two strings which ggplot2 recognizes as color names. These will set the colors of the increase and decrease bars of the plot. Default = c("violetred3", "turquoise3").
#' @param y_label A string for the y axis title. Default = "Source → Target".
#' @param x_label A string for the x axis title. Default = "Number".
#'
#' @returns A ggplot2 plot object.
visualize_differences <- function(data,
                                 title="",
                                 colors=c("violetred3", "turquoise3"),
                                 y_label="Source → Target",
                                 x_label="Number") {
    # Create the plot
    p <- ggplot(data, aes(x = number, y = pair_label, fill = number >= 0)) +
         geom_bar(stat = "identity") +
         scale_fill_manual(values = c("TRUE" = colors[1], "FALSE" = colors[2]),
                           labels = c("TRUE" = "Increase", "FALSE" = "Decrease"),
                           name = "Changes in Nr of CCC events \nbetween two celltypes") +
         theme_minimal() +
         theme(axis.text.y = element_text(hjust = 0, size = 10, margin=margin(t=100, b=100), face="bold", family = "mono"),
               panel.grid.major.x = element_blank(),
               plot.title = element_text(size = 14,
                                         # Center title
                                         hjust = 0.5,
                                         # Add space below title
                                         margin = margin(b = 10)),
               plot.title.position = "plot") +
        labs(title = title,
             y = y_label,
             x = x_label) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black")
    
    # Return the plot
    return(p)
}



#' Prepares data for the visualize_differences function. Two dataframes containing source celltype, target celltype and number of interactions between them for one condition/timepoint in a dataset are mergen to get the difference in interaction numbers between celltype pairs. 
#' Author: Krsitina Müller (kmlr81)
#'
#' @param data_frame1, data_frame2 Two data frames of the same structure. One per condition to be compared. 
#' @returns A data frame with the same three columns as the input data frames but now the number values show changes in condition 1 compared to condition 2. (df1 - df2).
prep_data_barplot <- function(data_frame1, data_frame2) {
    result <- merge(data_frame1, data_frame2, by = c("source_cluster", "target_cluster"), suffixes = c("_df1", "_df2"), all = TRUE)

    # Names for condition specific number columns
    number_col1 <- "number_df1"
    number_col2 <- "number_df2"

    # Subtract number_col2 from number_col1 
    # handle cases where source/target combination isn't present in one of the conditions
    result$number <- ifelse(!is.na(result[[number_col1]]) & !is.na(result[[number_col2]]), 
                            result[[number_col1]] - result[[number_col2]],
                            ifelse(!is.na(result[[number_col1]]), 
                                   result[[number_col1]],
                                   -result[[number_col2]])
                           )
    # Keep only relevant columns
    result <- result[, c("source_cluster", "target_cluster", "number")]

    # Create a new column for y labels
    max_length <- max(nchar(result$source))
    result$source_padded <- sapply(result$source_cluster, function(x) {paste0(x, strrep(" ", max_length - nchar(x)))})
    result$pair_label <- paste0(result$source_padded, "  →  ", result$target_cluster)
    #result <- result %>% mutate(pair_label = paste(source_cluster, "→", target_cluster)) %>%
              #mutate(pair_label = factor(pair_label, levels = pair_label[order(source_cluster, target_cluster, decreasing = TRUE)]))
    
    # Return new data frame
    return(result)
}


#' Runs SoupX background correction on given count matrix (filtered)
#' Author: Krsitina Müller (kmlr81) 
#' (code from chapter 6.4 of the book [Single-cell best practices](https://www.sc-best-practices.org/preamble.html) after Heumos, L.)
#'
#' @param cmatrix_filtered  A count matrix with empty cells filtered out
#' @param cmatrix_raw A count matrix with empty cells 
#' @param genes A list of gene names for the matrices
#' @param cells A list of cell tags for the matrices
#' @param soupx_groups A list of cluster labels from basic leiden clustering for the data
#' @param plot_path A string with the path to save the plot of the SoupX background correction
#' @param forceAccept A boolean to force the acceptance of the SoupX background correction. Default = FALSE
#' @returns cmatrix_soupx_out backgroundcorrected count matrix
run_SoupX <- function(cmatrix_filtered,
                      cmatrix_raw,
                      genes,
                      cells,
                      soupx_groups,
                      plot_path,
                      forceAccept=FALSE) {

    # Set row names and col names for the filtered count matrix
    rownames(cmatrix_filtered) = genes
    colnames(cmatrix_filtered) = cells

    # Set correct sparese format for count matrices
    cmatrix_filtered <- as(cmatrix_filtered, "sparseMatrix")
    cmatrix_raw <- as(cmatrix_raw, "sparseMatrix")

    # Generate SoupChannel object for SoupX
    soup_channel = SoupChannel(cmatrix_raw, cmatrix_filtered, calcSoupProfile=FALSE)
    # Add additional meta data to the SoupChannel object
    soup_profile = data.frame(row.names=rownames(cmatrix_filtered), est=rowSums(cmatrix_filtered)/sum(cmatrix_filtered), counts=rowSums(cmatrix_filtered))
    soup_channel = setSoupProfile(soup_channel, soup_profile)
    # Set cluster information for SoupChannel object
    soup_channel = setClusters(soup_channel, soupx_groups)

    png(plot_path)
    soup_channel = autoEstCont(soup_channel, forceAccept=forceAccept)
    dev.off()

    # Get corrected count matrix and round to counts Integer
    cmatrix_soupx_out = adjustCounts(soup_channel, roundToInt=TRUE)

    return(cmatrix_soupx_out)
}


#' Generate the annotations for the complex heatmap overview of ligand-target regulatory potential
#' with AUPR, AUROC and mean gene expression scores per sender cell type for the top ligands
#' Author: Kristina Müller, Example code that was modified for this function by: Holger Kristen
#' 
#' @param ligand_activities data frame with nichenetr scores for ligands
#' @param top_ligand string vecotr of the names of top ligands to consider when plotting
#' sender_celltypes string vector with names of all potential sender cell types
#' mean expressions matrix with mean expressions per sender cell type for each top ligand
#' @returns a row annotation object for the ComplexHeatmap Heatmap() function
generate_heatmap_annotations <- function(ligand_activities,
                                         top_ligands,
                                         sender_celltypes,
                                         mean_expressions) {

    # Get auroc data frame of top ranked ligands
    vis_ligand_auroc <- ligand_activities %>% filter(test_ligand %in% top_ligands) %>%
    column_to_rownames("test_ligand") %>% select(auroc) %>% as.matrix(ncol = 1)
    vis_ligand_auroc <- vis_ligand_auroc[top_ligands, ]

    # Get aupr data frame of top ranked ligands
    vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% top_ligands) %>%
    column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% as.matrix(ncol = 1)
    vis_ligand_aupr <- vis_ligand_aupr[top_ligands, ]

    # Replace unwanted characters in sender cell type
    bad_chars <- c("/", " ")
    clean_sender_celltypes <- gsub(paste(bad_chars, collapse = "|"), "-", sender_celltypes)

    # Prepare mean expression data frame for plotting
    # Noramlize each row of the mean expressiond data frame to values between 0 - 1 
    mean_expr_df_norm <- as.data.frame(t(apply(mean_expressions, 1, function(x) {(x - min(x)) / (max(x) - min(x))})))
    
    # Transpose and order ligands
    mean_expr_df_t_norm <- t(mean_expr_df_norm)[top_ligands,]
    colnames(mean_expr_df_t_norm) <- clean_sender_celltypes

    # Transpose and order ligands for non-normalized expression data frame 
    #mean_expr_df_t <- t(mean_expressions)[top_ligands,]
    #colnames(mean_expr_df_t) <- clean_sender_celltypes

    # Get max value for mean gene expression across all celltypes for color scale
    mean_max_norm <- max(mean_expr_df_t_norm)
    #mean_max <- max(mean_expr_df_t)
    
    # Create row annotation
    row_anno <- rowAnnotation(AUPR=row_anno_barplot(vis_ligand_aupr), 
                              "AUROC-0.5"=row_anno_barplot(vis_ligand_auroc - 0.5),
                              senders_norm=mean_expr_df_t_norm,
                              col=list(senders_norm=colorRamp2(c(0, round(mean_max_norm, digits=2)), c("white", "purple"))),
                              annotation_name_rot=75,
                              annotation_legend_param=list(senders_norm=list(title="Mean gene expression \nin sender celltpyes \nnormalized")))

    # return row annotations object
    return(row_anno)
}


#' Generates and saves a complexHeatmap object of ligand target interaction potential with 
#' row annotations showing the normalized mean expression of all ligands across sender cell types, aupr and auroc scores for ligands.
#' Function by: Kristina Müller, kmlr81
#'
#' @param matrix A matrix with the regulatory potential of ligand-target interactions
#' @param right_annotation A row annotation object for the ComplexHeatmap Heatmap() function
#' @param file_path A string with the path to save the heatmap plot to as a pdf. Default = "". 
#'                  If default the plot will be saved to current working directory.
#' @param file_width A number for the width of the pdf plot. Default = 13.
#' @param file_height A number for the height of the pdf plot. Default = 10.
#' @param row_order A vector with the order of rows in the heatmap. Default = NULL.
#' @param cluster_rows A boolean that decides whether to cluster rows in the heatmap. Default = FALSE.
#' @param column_order A vector with the order of columns in the heatmap. Default = NULL.
#' @param cluster_columns A boolean that decides whether to cluster columns in the heatmap. Default = FALSE.
#' @param row_dend_width A unit object for the width of the row dendrogram. Default = unit(20, "cm").
#' @param row_names_side A string for the side of the row names in the plot. Default = "left".
#' @param column_dend_height A unit object for the height of the column dendrogram. Default = unit(1.5, "cm").
#' @param column_title A string for the title of the columns in the plot. Default = "Target genes".
#' @param column_title_side A string for the side of the column title in the plot. Default = "bottom".
#' @param column_title_gp A gpar object for formatting options of the column title in the plot. Default = gpar(fontsize=15).
#' @param row_title A string for the title of the rows in the plot. Default = "Top ligands".
#' @param row_title_side A string for the side of the row title in the plot. Default = "left".
#' @param row_title_gp A gpar object for formatting options of the row title in the plot. Default = gpar(fontsize=15).
#' @param col A colorRamp2 object for the color scale of the heatmap. Default = colorRamp2(c(0, 1), c("white", "green")).
#' @param name A string for the name of the heatmap. Default = "complexHeatmap".
#' @param heatmap_legend_param A list with parameters for the heatmap legend. Default = list(title="Regulatory \npotential").
#' @param heatmap_padding A unit object for the padding of the heatmap. Default = unit(c(2,2,20,2), "mm").
#' @param heatmap_title A string for the title of the heatmap. Default = "".
#' @param heatmap_title_y A unit object for the y position of the heatmap title. Default = unit(1, "npc") + unit(2, "mm").
#' @param heatmap_title_just A string for the position of the heatmap title. Default = "bottom".
#' @param heatmap_title_gp A gpar object for formatting options of the heatmap title. Default = gpar(fontsize=15).
generate_complex_heatmap <- function(matrix,
                                     right_annotation,
                                     file_path="",
                                     file_width=13,
                                     file_height=10,
                                     row_order=NULL,
                                     cluster_rows=FALSE,
                                     column_order=NULL,
                                     cluster_columns=FALSE,
                                     row_dend_width=unit(20, "cm"),
                                     row_names_side="left",
                                     column_dend_height=unit(1.5, "cm"),
                                     column_title="Target genes",
                                     column_title_side="bottom",
                                     column_title_gp=gpar(fontsize=15),
                                     row_title="Top ligands",
                                     row_title_side="left",
                                     row_title_gp=gpar(fontsize=15),
                                     col=colorRamp2(c(0, 1), c("white", "green")),
                                     name="complexHeatmap",
                                     heatmap_legend_param=list(title="Regulatory \npotential"),
                                     heatmap_padding=unit(c(2,2,20,2), "mm"),
                                     heatmap_title="",
                                     heatmap_title_y=unit(1, "npc") + unit(2, "mm"),
                                     heatmap_title_just="bottom",
                                     heatmap_title_gp=gpar(fontsize=15)){
    pdf(file_path, width=file_width, height=file_height)


    heatmap <- Heatmap(matrix,
                       row_order=row_order,
                       cluster_rows=ifelse(is.null(row_order), TRUE, cluster_rows), 
                       row_dend_width=row_dend_width,
                       row_names_side=row_names_side,
                       column_order=column_order,
                       cluster_columns=ifelse(is.null(column_order), TRUE, cluster_columns),
                       column_dend_height=column_dend_height, 
                       column_title=column_title,
                       column_title_side=column_title_side,
                       column_title_gp=column_title_gp,
                       row_title=row_title,
                       row_title_side=row_title_side,
                       row_title_gp=row_title_gp,
                       col=col,
                       name=name,
                       heatmap_legend_param=heatmap_legend_param, 
                       right_annotation=right_annotation)

    draw(heatmap, padding=heatmap_padding) # Add space for title with padding

    decorate_heatmap_body(name, {
        grid.text(heatmap_title, y=heatmap_title_y, just=heatmap_title_just, gp=heatmap_title_gp)
    })

    dev.off()
}

#' Generates one barplot showing the sum of normalized magnitude of cell-cell-communication events for each 
#' sender and target cell type per species.
#' Author: Kristina Müller, kmlr81
#' 
#' @param data_senders A data frame with the columns 'source', 'species', 'sum_magnitude_per_ligand', 
#'                     'maxsum_magnitude_rel_source' and 'group'.
#' @param data_targets A data frame with the columns 'target', 'species', 'sum_magnitude_per_receptor', 
#'                     'maxsum_magnitude_rel_target' and 'group'.
#' @param path_out A string with the path to save the barplots to as a pdf. Default = "". 
#'                 If default the plot will be saved to current working directory.  
#' @param experiment_prefix A string to be added to the file name of the barplot. Default = "".
#' @param file_width A number for the width of the pdf plot. Default = 16.
#' @param file_height A number for the hight of the pdf plot. Default = 7.
generate_mean_magnitude_barplot <- function(data_senders,
                                            data_targets,
                                            path_out = "",
                                            experiment_prefix = "",
                                            file_width = 16,
                                            file_height = 7) {
  # Define color palette for plotting
  NogpaletteReihe <-  c("#CB769E", "#DE639A", "#A85C85", "#0081AF", 
                        "#4F6D7A", "#7C6A0A", "#368F8B", "#246A73", 
                        "#5CC1BC", "#62C370", "#F7C548", "#F97E44", 
                        "#FB3640", "#B7245C", "#0D3B66", "#3E2F5B", 
                        "#B2675E", "#644536")
  
  # Initialize empty lists for sender plots and target plots
  plots_senders <- list()
  plots_targets <- list()
  
  # Generate generate one senders and one targets plot per species with normalized magnitudes
  for (species in unique(data_senders$species)) {
    # Get senders and targets data frames for current species
    data_species_senders <- data_senders[data_senders$species == species, ]
    data_species_targets <- data_targets[data_targets$species == species, ]
    
    # Generate senders bar plot for current species
    p_species_senders <- ggplot(data_species_senders, 
                                aes(x=reorder(source %>% str_replace_all("_", " ") %>% str_wrap(width = 17), 
                                              -maxsum_magnitude_rel_source), 
                                    y=as.numeric(sum_magnitude_per_ligand), 
                                    fill=factor(group, levels=rev(unique(group))))) + 
      geom_col(position=position_dodge(width=0.9)) + 
      coord_flip() + 
      theme_minimal(base_size=16) + 
      xlab("") + 
      ylab("") + 
      scale_fill_manual(values=sample(NogpaletteReihe, length(unique(data_species_senders$group)), replace=FALSE), 
                        guide=guide_legend(reverse=TRUE, title=NULL)) + 
      ggtitle(species) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    
    # Generate targets bar plot for current species
    p_species_targets <- ggplot(data_species_targets, 
                                aes(x=reorder(target %>% str_replace_all("_", " ") %>% str_wrap(width = 17), 
                                              -maxsum_magnitude_rel_target), 
                                    y=as.numeric(sum_magnitude_per_receptor), 
                                    fill=factor(group, levels=rev(unique(group))))) + 
      geom_col(position=position_dodge(width=0.9)) + 
      coord_flip() + 
      theme_minimal(base_size=16) + 
      xlab("") + 
      ylab("") + 
      scale_fill_manual(values=sample(NogpaletteReihe, length(unique(data_species_targets$group)), replace=FALSE), 
                        guide=guide_legend(reverse=TRUE, title=NULL)) + 
      ggtitle(species) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    
    # Add plots to lists
    plots_senders <- c(plots_senders, list(p_species_senders))
    plots_targets <- c(plots_targets, list(p_species_targets))
  }
  
  # Ensure y-axis ticks are shown for the first plot in senders and targets plots
  plots_senders[[1]] <- plots_senders[[1]] + theme(axis.text.y = element_text(), axis.ticks.y = element_line())
  plots_targets[[1]] <- plots_targets[[1]] + theme(axis.text.y = element_text(), axis.ticks.y = element_line())
  
  # combine all senders plots for species into one plot and add caption
  p_all_senders <- wrap_plots(plots_senders) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = "Sum of normalized magnitude of cell-cell-communication events - \nsenders (ligands)"
    ) &
    theme(
      plot.caption = element_text(hjust = 0.5, size = 16),
      legend.position = "right"
    )

  # Combine all targets plots for speices into one plot and add caption
  p_all_targets <- wrap_plots(plots_targets) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = "Sum of normalized magnitude of cell-cell-communication events - \ntargets (receptors)"
    ) &
    theme(
      plot.caption = element_text(hjust = 0.5, size = 16),
      legend.position = "right"
    )

  # Save combined plots for senders and targets to pdf files
  pdf(file.path(path_out, paste0(experiment_prefix, "_magnitude_senders_grouped.pdf")), width = file_width, height = file_height)
  print(p_all_senders)
  dev.off()

  pdf(file.path(path_out, paste0(experiment_prefix, "_magnitude_targets_grouped.pdf")), width = file_width, height = file_height)
  print(p_all_targets)
  dev.off()
}