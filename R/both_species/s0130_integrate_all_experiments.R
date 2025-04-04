#' ---
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     code_folding: hide
#' ---

# # Initialize -----
knitr::opts_chunk$set(cache = F, fig.width = 16, fig.height = 16)

# Install necessary packages if they are not already installed
require(toolboxH)
require(here)
library(Seurat)
library(toolboxH)
library(patchwork)
library(ggplot2)

# # load ====
seurat_cyno = readRDS( here("analysisR/results/s0125_seurat_QCed_cyno_orthologyzed1to1.rds"))
seurat_human = readRDS( here("analysisR/results/s0120_seurat_QCed_human.rds"))
ortho3_1to1 <- fread(here("analysisR/results/s0120_cyno2human_1to1_orthologues.csv.gz"))

qlist1 = venn2(colnames(seurat_cyno), colnames(seurat_human))
stopifnot(length(qlist1$q1) == 0)

qlist2 = venn2(rownames(seurat_cyno), rownames(seurat_human))
qlist3 = venn3(rownames(seurat_cyno), rownames(seurat_human), ortho3_1to1$`Human gene name`)



seurat2 = merge(seurat_human[qlist2$q1,], seurat_cyno[qlist2$q1,]) # , add.cell.ids = c("H24", "M24"), project = "H24_M24" not necessary as no cell id overlap
seurat2
# seurat2 <- RenameAssays(seurat2, originalexp = "RNA")
seurat2

seurat2$run10x %>% mytable

seurat2@assays
seurat2 <- JoinLayers(seurat2)


seurat2[["RNA"]] <- split(seurat2[["RNA"]], f = seurat2$run10x)
seurat2

seurat2 <- NormalizeData(seurat2)
seurat2 <- FindVariableFeatures(seurat2)
seurat2 <- ScaleData(seurat2)
seurat2 <- RunPCA(seurat2)

seurat2 <- FindNeighbors(seurat2, dims = 1:30, reduction = "pca")


# ## Unintegrated merged UMAP----

seurat2 <- FindClusters(seurat2, resolution = 2, cluster.name = "unintegrated_clusters")

seurat2 <- RunUMAP(seurat2, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

p_unintegrated = (DimPlot(seurat2, group.by = "run10x", label = T, reduction = "umap.unintegrated")  +
                    (DimPlot(seurat2, group.by = "louvain_res_7", label = T, reduction = "umap.unintegrated") + NoLegend()) ) / 
  (DimPlot(seurat2, group.by = "predicted.celltype.l1", label = T, reduction = "umap.unintegrated" ) +
     DimPlot(seurat2, group.by = "seurat_clusters", label = T, reduction = "umap.unintegrated") + NoLegend()) + plot_annotation(title = "Unintegrated merged data")


p_unintegrated

# # speichern unintegrated ----
pdf(here("analysisR/results/s0130_UMAP_cyno_human_merged_unintegrated.pdf"), 16,16)
p_unintegrated
dev.off()

saveRDS(seurat2, here("analysisR/results/s0130_merged_seurat2_unintegrated.rds"))





# ## Integrated merged , 5 variants----
# https://satijalab.org/seurat/articles/seurat5_integration

options(future.globals.maxSize = 10e+09)
options('future.globals.maxSize')
require(SeuratWrappers) # for FastMNNIntegration and scvi integration

DefaultAssay(seurat2)


seurat2 <- IntegrateLayers(
  object = seurat2, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
seurat2 <- IntegrateLayers(
  object = seurat2, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
seurat2 <- IntegrateLayers(
  object = seurat2, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
) # not available, also not on github https://raw.githubusercontent.com/satijalab/seurat/1549dcb3075eaeac01c925c4b4bb73c73450fc50/R/integration5.R, also scVIIntegration not available any more
seurat2 <- IntegrateLayers(
  object = seurat2, method = JointPCAIntegration,
  new.reduction = "integrated.jpca",verbose = FALSE
)

seurat2 <- IntegrateLayers(
  object = seurat2, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "C:/Users/hol7525kir/AppData/Local/r-miniconda/envs/scvi-env", verbose = FALSE ) # from windows anaconda prompt, befehl  conda info --envs, dann conda activate C:/Users/hol7525kir/AppData/Local/r-miniconda, dann pip install scvi-tools, dann update SeuratWrappers  . On server disc not possible filename too long

seurat2 <- IntegrateLayers(
  object = seurat2, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)



# # visualisieren ----
seurat2 <- FindNeighbors(seurat2, reduction = "integrated.rpca", dims = 1:30)
seurat2 <- FindClusters(seurat2, resolution = 2, cluster.name = "rpca_clusters")
seurat2 <- RunUMAP(seurat2, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

seurat2 <- FindNeighbors(seurat2, reduction = "harmony", dims = 1:30)
seurat2 <- FindClusters(seurat2, resolution = 2, cluster.name = "harmony_clusters")
seurat2 <- RunUMAP(seurat2, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")


seurat2 <- FindNeighbors(seurat2, reduction = "integrated.mnn", dims = 1:30)
seurat2 <- FindClusters(seurat2, resolution = 2, cluster.name = "mnn_clusters")
seurat2 <- RunUMAP(seurat2, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")

seurat2 <- FindNeighbors(seurat2, reduction = "integrated.jpca", dims = 1:30)
seurat2 <- FindClusters(seurat2, resolution = 2, cluster.name = "jpca_clusters")
seurat2 <- RunUMAP(seurat2, reduction = "integrated.jpca", dims = 1:30, reduction.name = "umap.jpca")


seurat2 <- FindNeighbors(seurat2, reduction = "integrated.scvi", dims = 1:30)
seurat2 <- FindClusters(seurat2, resolution = 2, cluster.name = "scvi_clusters")
seurat2 <- RunUMAP(seurat2, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")


seurat2 <- FindNeighbors(seurat2, reduction = "integrated.cca", dims = 1:30)
seurat2 <- FindClusters(seurat2, resolution = 2, cluster.name = "cca_clusters")
seurat2 <- RunUMAP(seurat2, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")


p_integrated = DimPlot(seurat2, group.by = "celltypes_v0", label = T, reduction = "umap.rpca")  +
  DimPlot(seurat2, group.by = "celltypes_v0", label = T, reduction = "umap.harmony" ) +
  DimPlot(seurat2, group.by = "celltypes_v0", label = T, reduction = "umap.mnn") + 
  DimPlot(seurat2, group.by = "celltypes_v0", label = T, reduction = "umap.jpca") +
  DimPlot(seurat2, group.by = "celltypes_v0", label = T, reduction = "umap.scvi")  +
  DimPlot(seurat2, group.by = "celltypes_v0", label = T, reduction = "umap.cca")  + plot_annotation(title = "Integrated merged data") + plot_layout(ncol = 3, guides = "collect" ) & theme(legend.position = "right") 


p_integrated

p_integrated2 = DimPlot(seurat2, group.by = "rpca_clusters", label = T, reduction = "umap.rpca")  +
  DimPlot(seurat2, group.by = "harmony_clusters", label = T, reduction = "umap.harmony" ) +
  DimPlot(seurat2, group.by = "mnn_clusters", label = T, reduction = "umap.mnn") + 
  DimPlot(seurat2, group.by = "jpca_clusters", label = T, reduction = "umap.jpca") +
  DimPlot(seurat2, group.by = "scvi_clusters", label = T, reduction = "umap.scvi")  +
  DimPlot(seurat2, group.by = "cca_clusters", label = T, reduction = "umap.cca")  + plot_annotation(title = "Integrated merged data") + plot_layout(ncol = 3, guides = "collect" ) & guides(col = 'none') 


p_integrated2



p_integrated3 = DimPlot(seurat2, group.by = "predicted.celltype.l2", label = T, reduction = "umap.rpca")  +
  DimPlot(seurat2, group.by = "predicted.celltype.l2", label = T, reduction = "umap.harmony" ) +
  DimPlot(seurat2, group.by = "predicted.celltype.l2", label = T, reduction = "umap.mnn") + 
  DimPlot(seurat2, group.by = "predicted.celltype.l2", label = T, reduction = "umap.jpca") +
  DimPlot(seurat2, group.by = "predicted.celltype.l2", label = T, reduction = "umap.scvi")  +
  DimPlot(seurat2, group.by = "predicted.celltype.l2", label = T, reduction = "umap.cca")  + plot_annotation(title = "Integrated merged data") + plot_layout(ncol = 3, guides = "collect" ) & guides(col = 'none') 


p_integrated3

p_integrated4 =  (VlnPlot(seurat2,"IL2RA", group.by = "rpca_clusters", pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'rpca_clusters') )+
 ( VlnPlot(seurat2,"IL2RA",  group.by = "harmony_clusters" , pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'harmony_clusters')) +
  (VlnPlot(seurat2,"IL2RA",  group.by = "mnn_clusters" , pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'mnn_clusters') ) + 
  (VlnPlot(seurat2,"IL2RA",  group.by = "jpca_clusters" , pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'jpca_clusters') ) +
 ( VlnPlot(seurat2,"IL2RA",  group.by = "scvi_clusters" , pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'scvi_clusters'))  +
  (VlnPlot(seurat2,"IL2RA",  group.by = "cca_clusters", pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'cca_clusters')  ) + plot_annotation(title = "Integrated merged data") + plot_layout(ncol = 3, guides = "collect" ) & guides(col = 'none') 


activ_celltypes = grep("Prolif", unique(seurat2$predicted.celltype.l2), value = T)
activ_celltypes 
activ_celltypes = c("CD4 Proliferating", "CD8 Proliferating")
p_integrated4

seurat3 =  JoinLayers(seurat2)

seurat3activated = subset(seurat3,  predicted.celltype.l2 %in% activ_celltypes)


p_integrated5 =  (VlnPlot(seurat3activated,"IL2RA", group.by = "rpca_clusters", split.by =  "predicted.celltype.l2" ,pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'rpca_clusters') )+
  ( VlnPlot(seurat3activated,"IL2RA",  group.by = "harmony_clusters" ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'harmony_clusters')) +
  (VlnPlot(seurat3activated,"IL2RA",  group.by = "mnn_clusters" ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'mnn_clusters') ) + 
  (VlnPlot(seurat3activated,"IL2RA",  group.by = "jpca_clusters"  ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'jpca_clusters') ) +
  ( VlnPlot(seurat3activated,"IL2RA",  group.by = "scvi_clusters"  ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'scvi_clusters'))  +
  (VlnPlot(seurat3activated,"IL2RA",  group.by = "cca_clusters" ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("IL2RA", 'cca_clusters')  ) + plot_annotation(title = "Integrated merged data") + plot_layout(ncol = 3, guides = "collect" ) & guides(col = 'none') 


p_integrated5


p_integrated6 =  (VlnPlot(seurat3activated,"MKI67", group.by = "rpca_clusters", split.by =  "predicted.celltype.l2" ,pt.size = 0.1, alpha= 0.2) +ggtitle("MKI67", 'rpca_clusters') )+
  ( VlnPlot(seurat3activated,"MKI67",  group.by = "harmony_clusters" ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("MKI67", 'harmony_clusters')) +
  (VlnPlot(seurat3activated,"MKI67",  group.by = "mnn_clusters" ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("MKI67", 'mnn_clusters') ) + 
  (VlnPlot(seurat3activated,"MKI67",  group.by = "jpca_clusters"  ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("MKI67", 'jpca_clusters') ) +
  ( VlnPlot(seurat3activated,"MKI67",  group.by = "scvi_clusters"  ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("MKI67", 'scvi_clusters'))  +
  (VlnPlot(seurat3activated,"MKI67",  group.by = "cca_clusters" ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("MKI67", 'cca_clusters')  ) + plot_annotation(title = "Integrated merged data") + plot_layout(ncol = 3, guides = "collect" ) & guides(col = 'none') 


p_integrated6


p_integrated7 =  (VlnPlot(seurat3activated,"TFRC", group.by = "rpca_clusters", split.by =  "predicted.celltype.l2" ,pt.size = 0.1, alpha= 0.2) +ggtitle("TFRC", 'rpca_clusters') )+
  ( VlnPlot(seurat3activated,"TFRC",  group.by = "harmony_clusters" ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("TFRC", 'harmony_clusters')) +
  (VlnPlot(seurat3activated,"TFRC",  group.by = "mnn_clusters" ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("TFRC", 'mnn_clusters') ) + 
  (VlnPlot(seurat3activated,"TFRC",  group.by = "jpca_clusters"  ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("TFRC", 'jpca_clusters') ) +
  ( VlnPlot(seurat3activated,"TFRC",  group.by = "scvi_clusters"  ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("TFRC", 'scvi_clusters'))  +
  (VlnPlot(seurat3activated,"TFRC",  group.by = "cca_clusters" ,split.by =  "predicted.celltype.l2", pt.size = 0.1, alpha= 0.2) +ggtitle("TFRC", 'cca_clusters')  ) + plot_annotation(title = "Integrated merged data") + plot_layout(ncol = 3, guides = "collect" ) & guides(col = 'none') 


p_integrated7


p_integrated8 =  FeaturePlot(seurat2, "IL2RA", reduction = "umap.rpca")  +
  FeaturePlot(seurat2, "IL2RA", reduction = "umap.harmony" ) +
  FeaturePlot(seurat2, "IL2RA" ,reduction = "umap.mnn") + 
  FeaturePlot(seurat2, "IL2RA", reduction = "umap.jpca") +
  FeaturePlot(seurat2, "IL2RA", reduction = "umap.scvi")  +
  FeaturePlot(seurat2, "IL2RA" ,reduction = "umap.cca")  + plot_annotation(title = "Integrated merged data") + plot_layout(ncol = 3, guides = "collect" ) & guides(col = 'none') 

p_integrated8


p_integrated9 =  FeaturePlot(seurat2, "MKI67", reduction = "umap.rpca")  +
  FeaturePlot(seurat2, "MKI67", reduction = "umap.harmony" ) +
  FeaturePlot(seurat2, "MKI67" ,reduction = "umap.mnn") + 
  FeaturePlot(seurat2, "MKI67", reduction = "umap.jpca") +
  FeaturePlot(seurat2, "MKI67", reduction = "umap.scvi")  +
  FeaturePlot(seurat2, "MKI67" ,reduction = "umap.cca")  + plot_annotation(title = "Integrated merged data") + plot_layout(ncol = 3, guides = "collect" ) & guides(col = 'none') 

p_integrated9

pdf(here("analysisR/results/s0130_first_azimuth_six_integration.pdf"), height = 16,23)
p_integrated
p_integrated2
p_integrated3
p_integrated4
p_integrated5
p_integrated6
p_integrated7
p_integrated8
p_integrated9
dev.off()

# copy pdf to results_GIT
file.copy(from = here("analysisR/results/s0130_first_azimuth_six_integration.pdf"), to = here("analysisR/results_GIT/s0130_first_azimuth_six_integration.pdf"), overwrite = TRUE)

# # speichern ----

saveRDS(seurat2, here("analysisR/results/s0130_merged_seurat2_integrated.rds"))


# # finalize ----
finalizeSkript()
