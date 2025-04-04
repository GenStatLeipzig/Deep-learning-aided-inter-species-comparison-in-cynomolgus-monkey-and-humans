library(Seurat)
library(data.table)
library(ggplot2)
library(here)
library(stringr)

# Define marker groups
marker_groups_peter <- rbind(
  data.table(celltype = "B cells", markers = c("CD79B")),
  data.table(celltype = "Plasma cells", markers = c("PRDM1")),
  data.table(celltype = "T cells", markers = c("CD3E", "CD4", "CD8A")),
  data.table(celltype = "gdT", markers = c("GATA3")),
  data.table(celltype = "MAIT", markers = c("SLC4A10", "RORC")),
  data.table(celltype = "NK cells", markers = c("NCAM1", "KIR3DL1", "KLRD1")),
  data.table(celltype = "Mono", markers = c("FCN1", "CD14", "S100A8")),
  data.table(celltype = "Mono-Nonclass", markers = c("CX3CR1", "FCGR3A")),
  data.table(celltype = "DC", markers = c("FLT3", "FCER1A")),
  data.table(celltype = "Platelets", markers = c("PPBP", "PF4"))#,
  # data.table(celltype = "Proliferating", markers = c("TOP2A", "MKI67")),
  # data.table(celltype = "Neutrophils", markers = c("CXCR2"))
)

# Define the dot plot function
doMarkerDotPlot <- function(seurat, marker_groups_peter, grouping_factor) {
  marker_groups2 = copy(marker_groups_peter)
  marker_groups2[, highmarker := markers]
  marker_groups2[, celltype2 := celltype]
  marker_groups = copy(marker_groups2)
  
  marker_groups2[, highmarker := stringr::str_trim(highmarker)]
  
  # Create the base dotplot
  plot_data <- DotPlot(
    seurat,
    group.by = grouping_factor,
    features = unique(marker_groups2$highmarker),
    scale.by = "radius",
    scale = TRUE
  )
  
  # Extract the plot data
  plot_df <- as.data.table(plot_data$data, keep.rownames = TRUE)
  
  # Add grouping information
  plot_df2 = merge(
    plot_df, 
    marker_groups2[, .(highmarker, group = celltype)] %>% unique(), 
    by.x = 'features.plot', 
    by.y = 'highmarker', 
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  
  plot_df2[, group := factor(group, levels = unique(marker_groups$celltype))]
  
  # Create the modified plot
  markerplot = ggplot(plot_df2, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0),
      panel.spacing = unit(1, "lines")
    ) +
    facet_grid(. ~ group, scales = "free", space = "free") +
    labs(
      x = "Cell Type", 
      y = "Genes",
      color = "Average Expression", 
      size = "Percent Expressed"
    ) +
    coord_flip()
  
  return(markerplot)
}

doMarkerDotPlot <- function(seurat, marker_groups_peter,grouping_factor =  "predicted.celltype.l1_1st") {
  # seurat = seurat_human2
  # grouping_factor =  "predicted.celltype.l1.5"
  marker_groups2 = copy(marker_groups_peter)
  marker_groups2[, highmarker := markers]
  marker_groups2[, celltype2 := celltype]
  marker_groups = copy(marker_groups2)
  
  
  marker_groups2[, highmarker := str_trim(highmarker)]
  marker_groups2[,.N, celltype]
  marker_groups2[,.N, celltype2] %>% data.frame()
  
  # qlist64a  = venn2(rownames(seurat_human), marker_groups2$highmarker)
  # qlist64b = venn2(rownames(seurat_cyno), marker_groups2$highmarker)
  
  
  # Create the base dotplot
  plot_data <- DotPlot(
    seurat,
    group.by = grouping_factor,
    # seurat_human2,
    #  group.by ="predicted.celltype.l2_1st",
    features = unique(marker_groups2$highmarker %>% unique()),
    scale.by = "radius",
    scale = TRUE)
  
  # Extract the plot data
  plot_df <- plot_data$data %>% as.data.table(keep.rownames = T)
  plot_df[,.N, id]
  
  # Add grouping information
  marker_groups2[, .(highmarker, celltype2)] %>% unique() %>% .[ toolboxH::allDuplicatedEntries(highmarker)]
  plot_df2 = merge(plot_df, marker_groups2[, .(highmarker, group = celltype)] %>% unique(), by.x = 'features.plot', by.y = 'highmarker', all.x = T,allow.cartesian=TRUE )
  
  
  plot_df2[,group := factor(group, levels = unique(marker_groups$celltype))]
  
  unique(plot_df2$id) %>% dput()
  
  levels_l1 = c("B", 
                "CD4 T", "CD8 T","other T", "NK",  "Mono", "DC", "other")
  levels_l2 =  c("B intermediate", "B memory", "B naive","Plasmablast", "CD4 Naive", 
                 "CD4 Proliferating", "CD4 TCM", "CD4 TEM","CD4 CTL" , "CD8 Naive",  "CD8 Proliferating" , "CD8 TCM", 
                 "CD8 TEM", "Treg",  "dnT", "gdT","MAIT", "ILC",  "NK", "NK Proliferating", "NK_CD56bright"  , "CD14 Mono", "CD16 Mono", "HSPC",  "Platelet", "cDC1", 
                 "cDC2", "pDC")
  
  levels_l1_5 = c("B", "Plasmablast", 
                  "CD4 T",  "CD4 Proliferating", "CD8 T",
                  "CD8 Proliferating", "dnT", "gdT", "MAIT", "other T", "NK", "NK Proliferating", "Mono", "CD14 Mono", "CD16 Mono", "DC","Platelet" ,"other")
  
  
  
  if(all(plot_df2$id %in%  levels_l2))  plot_df2[,id := factor(id, levels =  levels_l2)] else if(all(plot_df2$id %in%  levels_l1))  plot_df2[,id := factor(id, levels =  levels_l1)] else if(all(plot_df2$id %in%  levels_l1_5))  plot_df2[,id := factor(id, levels =  levels_l1_5)] else  if(all(plot_df2$id %in%  levels_l2))  plot_df2[,id := factor(id, levels =  c(levels_l1_5, levels_l2, levels_l1, plot_df2$id) %>% unique())]
  
  
  
  # venn2(levels_l1, levels_l2)
  
  # Create the modified plot
  markerplot = ggplot(plot_df2, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0),  # Left-align strip text
      # strip.placement = "outside",                          # Place strips outside
      panel.spacing = unit(1, "lines")                   # Adjust space between facets
    )+
    facet_grid(. ~ group, switch = "x",scales = "free", space = "free") +
    labs(x = "Cell Type", y = "Genes",
         color = "Average Expression", size = "Percent Expressed") +
    coord_flip() + ggtitle(deparse(substitute(seurat)), subtitle = deparse(substitute(grouping_factor)))
  
  print(markerplot)
  markerplot
}

# Load the Seurat object
seurat2 <-  readRDS(here("analysisR/results/s0140_seurat2.rds"))

seurat2$species |> table()

# Create species-specific plots
human_plot <- doMarkerDotPlot(
  seurat2[, seurat2$species == "human"], 
  marker_groups_peter, 
  "cluster_azimut1_5_scanvi_nkuni"
) + ggtitle("H.sapiens", subtitle = "")

human_plot

cyno_plot <- doMarkerDotPlot(
  seurat2[, seurat2$species == "cyno"], 
  marker_groups_peter, 
  "cluster_azimut1_5_scanvi_nkuni"
) + ggtitle("M.fascicularis", subtitle = "")
cyno_plot
# Combine plots
combined_plot <- cowplot::plot_grid(
  human_plot, cyno_plot,
  ncol = 1,
  labels = c("A", "B"),
  label_size = 16
)
combined_plot
# Save the plot
ggsave(
  here("analysisR/results_GIT//s0150_marker_dotplot_combined.pdf"),
  combined_plot,
  width = 8,
  height = 9,
  units = "in"
)
ggsave(
  here("analysisR/results//s0150_marker_dotplot_combined.pdf"),
  combined_plot,
  width = 8,
  height = 9,
  units = "in"
)
