require(Seurat)
require(here)
library(patchwork)
library(toolboxH)

ggplotSankey <- function(input, sort_by_frequency =T) {
  require(ggsankey)
  require(ggrepel)
  
  require(dplyr)
  sankey_gg_data = ggsankey::make_long(input, names(input))
  
  
  sankey_gg_data_dt = as.data.table(sankey_gg_data, keep.rownames = F)
  ordered = sankey_gg_data_dt[,.N, .(node, next_node)][order(is.na(next_node),-N, na.last = T)  ]
  # ordered2 = ordered[order(-node, -N)]
  
  # sankey_gg_data$next_node = factor(sankey_gg_data$next_node, levels = reihenfolge)
  NogpaletteReihe <-  c("#CB769E", "#DE639A", "#A85C85", "#0081AF", "#4F6D7A", "#7C6A0A", "#368F8B", "#246A73", "#5CC1BC", "#62C370", "#F7C548", "#F97E44", "#FB3640", "#B7245C", "#0D3B66", "#3E2F5B", "#B2675E", "#644536") %>% rev()
  
  
  if(sort_by_frequency ==T){
    sankey_gg_data$next_node = factor(sankey_gg_data$next_node, levels = unique(ordered$next_node) )
    sankey_gg_data$node = factor(sankey_gg_data$node, levels = unique(ordered$node))
  }
  p_sankey = ggplot(sankey_gg_data, aes(x = x, next_x = next_x,
                                        node = node, next_node = next_node,
                                        fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = alpha("grey55", 0.2), alpha = 0.4) +
    geom_sankey_text(size = 4, color = "black") +
    # scale_fill_manual(values = rep(NogpaletteReihe, 10) %>% sample()) +
    theme_sankey(base_size = 18) +
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5),
          axis.text.x = element_text(color = "black"))
  
  
  
  p_sankey
}


# Load the Seurat object
seurat2 <-  readRDS(here("analysisR/results/s0140_seurat2.rds"))

## add NK and NK prolif unified
seurat2$cluster_azimut1_5_scanvi_nkuni = fifelse(grepl("^NK",seurat2$cluster_azimut1_5_scanvi ), "NK", seurat2$cluster_azimut1_5_scanvi)
ggplotSankey(seurat2@meta.data %>% data.table %>% .[,.(cluster_azimut1_5_scanvi, cluster_azimut1_5_scanvi_nkuni)])



# Extract UMAP coordinates from integrated.scvi
umap_coords <- Embeddings(seurat2, reduction = "umap.scvi")
umap_dt <- as.data.table(umap_coords, keep.rownames=TRUE)
colnames(umap_dt) <- c('rn',"UMAP1", "UMAP2")

# Add cell type information
stopifnot(identical(umap_dt$rn, colnames(seurat2)))
umap_dt$predicted.celltype.l1.5 <- seurat2$predicted.celltype.l1.5
umap_dt$species <- seurat2$species
umap_dt$timepoint <- seurat2$timepoint
umap_dt$individual <- seurat2$individual
umap_dt$cluster_azimut1_5_scanvi_plusprolif = seurat2$cluster_azimut1_5_scanvi_plusprolif
umap_dt$cluster_azimut1_5_scanvi = seurat2$cluster_azimut1_5_scanvi 
umap_dt$cluster_azimut1_5_scanvi_nkuni = seurat2$cluster_azimut1_5_scanvi_nkuni 
umap_dt$scvi_clusters  = seurat2$scvi_clusters
# Create ggplot2 version vs. Seurat version
library(ggplot2)
p1 <- ggplot(umap_dt, aes(x = UMAP1, y = UMAP2, color = predicted.celltype.l1.5)) +
  geom_point(size = 0.5, alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  labs(title = "UMAP (ggplot2)",
       color = "Cell Type") +
  scale_color_discrete(guide = guide_legend(override.aes = list(size = 4)))

# Create Seurat version
p2 <- DimPlot(seurat2, 
              reduction = "umap.scvi",
              group.by = "predicted.celltype.l1.5",
              pt.size = 0.5,
              label = TRUE) +
  labs(title = "UMAP (Seurat)") +
  theme(legend.text = element_text(size = 8))

p1 + p2

# celltype different versions ----

p2b <- DimPlot(seurat2, 
              reduction = "umap.scvi",
              group.by = "cluster_azimut1_5_scanvi_plusprolif",
              pt.size = 0.5,
              label = TRUE) +
  labs(title = "UMAP (Seurat)") +
  theme(legend.text = element_text(size = 8)) + ggtitle("cluster_azimut1_5_scanvi_plusprolif")

p2c <-  DimPlot(seurat2, 
                  reduction = "umap.scvi",
                  group.by = "cluster_azimut1_5_scanvi",
                  pt.size = 0.5,
                  label = TRUE) +
  labs(title = "UMAP (Seurat)") +
  theme(legend.text = element_text(size = 8)) + ggtitle("cluster_azimut1_5_scanvi")

p2d <-  DimPlot(seurat2, 
                reduction = "umap.scvi",
                group.by = "cluster_azimut1_5_scanvi_nkuni",
                pt.size = 0.5,
                label = TRUE) +
  labs(title = "UMAP (Seurat)") +
  theme(legend.text = element_text(size = 8)) + ggtitle("cluster_azimut1_5_scanvi_nkuni")


(p2 +  ggtitle("predicted.celltype.l1.5") | p2b) / (p2c | p2d)

## plot for paper ----

pa = DimPlot(seurat2, 
        reduction = "umap.scvi",
        group.by = c("cluster_azimut1_5_scanvi_nkuni",'scvi_clusters' ),
        pt.size = 0.5,
        repel = TRUE,
        label = TRUE,label.size = 3) +
  labs(title = "UMAP (Seurat)") +
  theme(legend.text = element_text(size = 8)) & NoLegend()

pb = DimPlot(seurat2, 
             reduction = "umap.scvi",
             group.by = c('species', "timepoint"),
             pt.size = 0.5,
             label = FALSE) +
  labs(title = "UMAP (Seurat)") &
  theme(#legend.text = element_text(size = 8),
        legend.position = "top")  



# Display plots side by side




p_overview = pb + pa + plot_layout(widths = c(1,1,2.4)) & ggtitle("")
p_overview

ggsave(
  here("analysisR/results_GIT//s0160_Dimplot_UMAP_integrated.pdf"),
  p_overview,
  width = 14.5,
  height = 4,
  units = "in"
)
## subcluster CD8
  seurat_cd8 <-  seurat2[, seurat2$cluster_azimut1_5_scanvi == "CD8 T"]

pcd8_1=  DimPlot(seurat_cd8[,seurat_cd8$species=="cyno"], 
          reduction = "umap.scvi",
          group.by = c( 'scvi_clusters' ),
          split.by = "timepoint",
          ncol = 3,
          pt.size = 0.5,
          repel = TRUE,
          label = TRUE,label.size = 5) +
    labs(title = "") +
    theme(legend.text = element_text(size = 8)) & NoLegend()
  
pcd8_2=DimPlot(seurat_cd8[,seurat_cd8$species=="human"], 
        reduction = "umap.scvi",
        group.by = c( 'scvi_clusters' ),
        split.by = "timepoint",
        ncol = 3,
        pt.size = 0.5,
        repel = TRUE,
        label = TRUE,label.size = 5) +
  labs(title = "") +
  theme(legend.text = element_text(size = 8)) & NoLegend()

pcd8_1/pcd8_2 +plot_annotation(title = "CD8 T cell subclusters") 

pcd8_fin = pcd8_1/pcd8_2 +plot_annotation(title = "CD8 T cell subclusters") & ylim(-1,6) & theme_void() &  theme(strip.text.x = element_text(size = 12)) & NoLegend()
pcd8_fin

ggsave(
  here("analysisR/results_GIT//s0160_subclusterplot_CD8.pdf"),
  pcd8_fin,
  width = 6,
  height = 4,
  units = "in"
)

## subcluster CD16
sort(unique(seurat2$cluster_azimut1_5_scanvi))
seurat_cd16mono_human <-  subset(seurat2, subset = cluster_azimut1_5_scanvi == "CD16 Mono" & species == "human")
seurat_cd16mono_cyno <-  subset(seurat2, subset = cluster_azimut1_5_scanvi == "CD16 Mono" & species == "cyno")

pcd16mono_1=  DimPlot(seurat_cd16mono_cyno, 
                 reduction = "umap.scvi",
                 group.by = c( 'scvi_clusters' ),
                 split.by = "timepoint",
                 ncol = 3,
                 pt.size = 0.5,
                 repel = TRUE,
                 label = TRUE,label.size = 5) +
  labs(title = "") +
  theme(legend.text = element_text(size = 8)) & NoLegend()

timepoints_human <- factor(seurat_cd16mono_human$timepoint, 
                           levels = c("00hr", "06hr", "24hr"))
seurat_cd16mono_human$timepoint <- timepoints_human



seurat_cd16mono_human$timepoint %>% str()
pcd16mono_2=DimPlot(seurat_cd16mono_human, 
               reduction = "umap.scvi",
               group.by = c( 'scvi_clusters' ),
               split.by = "timepoint",
               ncol = 3,
                
               pt.size = 0.5,
               repel = TRUE,
               label = TRUE,label.size = 5) +
  labs(title = "") +
  theme(legend.text = element_text(size = 8)) +
  facet_wrap(~timepoint, ncol=3, drop=FALSE) &  NoLegend()

pcd16mono_1/pcd16mono_2 +plot_annotation(title = "cd16mono  subclusters")   &  theme(strip.text.x = element_text(size = 12)) & NoLegend()
pcd16mono_1/pcd16mono_2 +plot_annotation(title = "cd16mono  subclusters")  & theme_void() &  theme(strip.text.x = element_text(size = 12)) & NoLegend() & ylim(-8.5,-5.5)




## subcluster NK
sort(unique(seurat2$cluster_azimut1_5_scanvi))
seurat_NK_human <-  subset(seurat2, subset = cluster_azimut1_5_scanvi %in% c("NK", "NK Proliferating") & species == "human")
seurat_NK_cyno <-  subset(seurat2, subset = cluster_azimut1_5_scanvi  %in% c("NK", "NK Proliferating") & species == "cyno")

pNK_1=  DimPlot(seurat_NK_cyno, 
                      reduction = "umap.scvi",
                      group.by = c( 'scvi_clusters' ),
                      split.by = "timepoint",
                      ncol = 3,
                      pt.size = 0.5,
                      repel = TRUE,
                      label = TRUE,label.size = 5) +
  labs(title = "") +
  theme(legend.text = element_text(size = 8)) & NoLegend()

timepoints_human <- factor(seurat_NK_human$timepoint, 
                           levels = c("00hr", "06hr", "24hr"))
seurat_NK_human$timepoint <- timepoints_human



seurat_NK_human$timepoint %>% str()
pNK_2=DimPlot(seurat_NK_human, 
                    reduction = "umap.scvi",
                    group.by = c( 'scvi_clusters' ),
                    split.by = "timepoint",
                    ncol = 3,
                    
                    pt.size = 0.5,
                    repel = TRUE,
                    label = TRUE,label.size = 5) +
  labs(title = "") +
  theme(legend.text = element_text(size = 8)) +
  facet_wrap(~timepoint, ncol=3, drop=FALSE) &  NoLegend()

pNK_1/pNK_2 +plot_annotation(title = "NK  subclusters")   &  theme(strip.text.x = element_text(size = 12)) & NoLegend()
pNK_fin = pNK_1/pNK_2 +plot_annotation(title = "NK  subclusters")  & theme_void() &  theme(strip.text.x = element_text(size = 12)) & NoLegend() & ylim(-.5, 5)
pNK_fin

ggsave(
  here("analysisR/results_GIT//s0160_subclusterplot_NK.pdf"),
  pNK_fin,
  width = 6,
  height = 4,
  units = "in"
)

## speichern
fwrite(umap_dt, here("analysisR/results_GIT/s0160_UMAP_harmonized_scvi.txt.gz"))

## finalise
toolboxH::finalizeSkript()
