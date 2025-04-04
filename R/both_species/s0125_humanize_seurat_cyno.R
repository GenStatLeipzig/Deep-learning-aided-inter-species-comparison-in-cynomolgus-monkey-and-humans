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
# # renaming functions from verteszy github ----
require(Stringendo)

#seurat_ape to human gene names



# # Laden ----
# Load the Seurat object
seurat_cyno <- readRDS(here("analysisR/results/s0120_seurat_QCed_cyno.rds"))
seurat_cyno

ortho3_1to1 <- fread(here("analysisR/results/s0120_cyno2human_1to1_orthologues.csv.gz"))
ortho3_1to1[ `Gene name` != `Human gene name`]


showNA(ortho3_1to1) # ok, some have orthologue, but no human gene expression in my data and some same name, but no orthologue


renamer = ortho3_1to1[, .(seurat_gene, `Human gene name`)]
stopifnot(nrow(renamer[allDuplicatedEntries(seurat_gene)])==0)
stopifnot(nrow(renamer[allDuplicatedEntries(`Human gene name`)])==0)
renamer[seurat_gene=="PTPRC"]


renamer[seurat_gene != `Human gene name`]

# convert a v3 assay to a v5 assay
seurat_cyno2 = DietSeurat(seurat_cyno[renamer$seurat_gene,])
seurat_cyno2

countmatrix = seurat_cyno2@assays$RNA@counts 
hh(countmatrix,22)
dim(countmatrix)
str(countmatrix)


message("Renaming to orthologues :")

rownames(countmatrix) = renamer[match_hk(rownames(countmatrix), renamer$seurat_gene), `Human gene name`]

# Setup Seurat object - no filter
seurat_cyno2hum= CreateSeuratObject(countmatrix, min.cells = 1, min.features = 1) #  more liberal than standard as i will filter min.features, i.e. cells expressen a certain a certain number of features later
seurat_cyno2hum$orig.ident %>% table()
dim_after = dim(seurat_cyno2hum) %>% huebsch()

message("Importing seurat for with dimensions " , paste0(dim_after, collapse = " x "), " genes x cells.")  
seurat_cyno2hum


vars2add = setdiff( names(seurat_cyno@meta.data), names(seurat_cyno2hum@meta.data))
vars2add
# DefaultAssay(seurat_cyno2hum)<-"RNA"

for(i in vars2add){
  # i = vars2add[1]
  seurat_cyno2hum@meta.data[[i]] = seurat_cyno@meta.data[match_hk(rownames(seurat_cyno2hum@meta.data), rownames(seurat_cyno@meta.data)), i]  
  }

seurat_cyno2hum


qlist2 = venn3(ortho3_1to1$seurat_gene, rownames(seurat_cyno),rownames(seurat_cyno2hum))


(rownames(seurat_cyno2) == rownames(seurat_cyno2hum)) %>% mytable()

# # compare umap befor after humanization -----


seurat_cyno2hum <- NormalizeData(seurat_cyno2hum)
seurat_cyno2hum <- FindVariableFeatures(seurat_cyno2hum)
seurat_cyno2hum <- ScaleData(seurat_cyno2hum)
seurat_cyno2hum <- RunPCA(seurat_cyno2hum)

seurat_cyno2hum <- FindNeighbors(seurat_cyno2hum, dims = 1:30, reduction = "pca")
seurat_cyno2hum <- FindClusters(seurat_cyno2hum, resolution = 2, cluster.name = "unintegrated_clusters")

seurat_cyno2hum <- RunUMAP(seurat_cyno2hum, dims = 1:30)#, reduction = "pca", reduction.name = "umap.unintegrated")


seurat_cyno3 <- NormalizeData(DietSeurat(seurat_cyno))
# seurat_cyno3 <- NormalizeData(seurat_cyno2)
seurat_cyno3 <- FindVariableFeatures(seurat_cyno3)
seurat_cyno3 <- ScaleData(seurat_cyno3)
seurat_cyno3 <- RunPCA(seurat_cyno3)

seurat_cyno3 <- FindNeighbors(seurat_cyno3, dims = 1:30, reduction = "pca")
seurat_cyno3 <- FindClusters(seurat_cyno3, resolution = 2, cluster.name = "unintegrated_clusters")

seurat_cyno3 <- RunUMAP(seurat_cyno3, dims = 1:30)#, reduction = "pca", reduction.name = "umap.unintegrated")





p_compare = (DimPlot(seurat_cyno2hum, group.by = "run10x", label = T )  +
    (DimPlot(seurat_cyno2hum, group.by = "predicted.celltype.l1", label = T ) +
       VlnPlot(seurat_cyno2hum, "PTPRC", group.by =  "predicted.celltype.l1" )) ) / 
  (DimPlot(seurat_cyno3, group.by = "run10x", label = T )  +
     (DimPlot(seurat_cyno3, group.by = "predicted.celltype.l1", label = T ) +VlnPlot(seurat_cyno3, "PTPRC", group.by = "predicted.celltype.l1" )
      ))+ plot_annotation(title = "Original above vs. 1:1 Orthologue-only below")

p_compare






pdf(here("analysisR/results/s0125_UMAP_cyno_before_after_orthologues.pdf"), 16,16)
p_compare
dev.off()


pdf(here("analysisR/results_GIT/s0125_UMAP_cyno_before_after_orthologues.pdf"), 16,16)
p_compare
dev.off()

# # save ----

saveRDS(seurat_cyno2hum, here("analysisR/results/s0125_seurat_QCed_cyno_orthologyzed1to1.rds"))
saveRDS(seurat_cyno3, here("analysisR/results/s0125_seurat_QCed_cyno_NONorthologyzed.rds"))

# # finalize ----
finalizeSkript()