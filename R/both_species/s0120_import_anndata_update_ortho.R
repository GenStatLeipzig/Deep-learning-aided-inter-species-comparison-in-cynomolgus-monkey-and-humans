#' ---
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     code_folding: hide
#' ---



# # Initialize -----
knitr::opts_chunk$set(cache = F, fig.width = 16, fig.height = 12)

# Install necessary packages if they are not already installed
require(toolboxH)
require(here)
# Install zellkonverter and Seurat packages
# BiocManager::install("zellkonverter", lib = "C:/Users/hol7525kir/AppData/Local/R/win-library/4.2", force = T)  # ‘1.8.0
# remotes::install_github("theislab/zellkonverter") # ‘1.15.4’
# Load the required libraries
library(zellkonverter)
library(Seurat)
library(SummarizedExperiment)
library(SeuratObject)
library(toolboxH)
library(patchwork)

packageVersion("zellkonverter")

# # Read the AnnData as SingleCellExperiment object ----
sce_hum <- readH5AD(here("analysis/results/for_integration/H24_human_anno_celltypes_v0.h5ad"))
sce_hum
dim(sce_hum)
sce_ape <- readH5AD(here("analysis/results/for_integration/M24_cyno_anno_celltypes_v0.h5ad"))
sce_ape


# # Make seurat and visualise ----
# Check existing assays and keep X only, the most recent in python
SummarizedExperiment::assayNames(sce_hum)
sce_hum@assays@data$counts = NULL
sce_hum@assays@data$soupX_counts = NULL
assays(sce_hum)[["counts"]] <- assays(sce_hum)[["X"]]
sce_hum@assays@data$X = NULL

SummarizedExperiment::assayNames(sce_hum)
assayNames(sce_hum)

sce_hum2 =scuttle::logNormCounts(sce_hum)

seurat_hum <- as.Seurat(sce_hum2)
seurat_hum <- RenameAssays(seurat_hum, originalexp = "RNA")

(DimPlot(seurat_hum, group.by = "run10x", label = T)  + (DimPlot(seurat_hum, group.by = "louvain_res_7", label = T) + NoLegend()) ) / ((DimPlot(seurat_hum, group.by = "predicted.celltype.l1", label = T) ) + DimPlot(seurat_hum, group.by = "basic_QC", label = T) ) + plot_annotation(title = "H24_human_anno_celltypes_v0", subtitle = "LogNorm-based data layer added")

assayNames(sce_ape)
sce_ape@assays@data$counts = NULL
sce_ape@assays@data$soupX_counts = NULL
assays(sce_ape)[["counts"]] <- assays(sce_ape)[["X"]]
sce_ape@assays@data$X = NULL



assayNames(sce_ape)


assayNames(sce_ape)
assayNames(sce_ape)

sce_ape2 =scuttle::logNormCounts(sce_ape)

seurat_ape <- as.Seurat(sce_ape2)
seurat_ape <- RenameAssays(seurat_ape, originalexp = "RNA")

(DimPlot(seurat_ape, group.by = "run10x", label = T)  + (DimPlot(seurat_ape, group.by = "louvain_res_7", label = T) + NoLegend()) ) / ((DimPlot(seurat_ape, group.by = "predicted.celltype.l1", label = T) ) + DimPlot(seurat_ape, group.by = "basic_QC", label = T) ) + plot_annotation(title = "M24_cyno_anno_celltypes_v0", subtitle = "LogNorm-based data layer added")


# # Make a table including all seurat cyno genes annotated with orthologues ----
ortho = fread(here("preprocessing/ensembl/mart_export_orthologues.txt.gz"), na.strings = c("NA", "NULL", ""))
qlist1 = venn2(ortho$`Human gene name`, rownames(seurat_hum))

cynoanno= fread(here("preprocessing/ensembl/mart_export.txt.gz"), na.strings = c("NA", "NULL", ""))
qlist2 = venn3(cynoanno$`Gene name`,cynoanno$`Gene stable ID`, rownames(seurat_ape))
str(qlist2)
sort(qlist2$q7, decreasing = T)

qlist2b = venn3(cynoanno$`Gene name` %>% toupper(),cynoanno$`Gene stable ID`, rownames(seurat_ape) %>% toupper())
sort(qlist2b$q7, decreasing = T)

qlist3 = venn2(ortho$`Gene stable ID`, cynoanno$`Gene stable ID`)
qlist4 = venn2(names(ortho), names(cynoanno))

setkeyv(ortho, qlist4$q1)
setkeyv(cynoanno, qlist4$q1)


# sometimes multiple entry in cynoanno reagardind Ensembl Canonical, e.g. ZYG11A and sometimes multiple EntrezGene transcript name ID, e.g. ENSMFAG00000000229 and Gene Synonym, e.g. ENSMFAG00000046325
# collapsing this info

stopifnot(nrow(cynoanno[, -c('Ensembl Canonical', "EntrezGene transcript name ID", "Gene Synonym")] %>% unique() %>% .[allDuplicatedEntries(`Gene stable ID`)])==0)

cynoanno2 = cynoanno[, .(`Ensembl Canonical2` = paste(unique(sort(`Ensembl Canonical`)), collapse = ","),
                         `EntrezGene transcript name ID` = paste(unique(sort(`EntrezGene transcript name ID`)), collapse = ","),
                         `Gene Synonym2` = paste(unique(sort(`Gene Synonym`)), collapse = ",")
                         ), by = .(`Gene stable ID`,`Gene stable ID version`,`Gene start (bp)`,`Gene end (bp)`,`Chromosome/scaffold name`,`Gene description`,`Gene name`,`Source of gene name`,`Source (gene)`,`HGNC symbol`)]


stopifnot(nrow(cynoanno2[allDuplicatedEntries(`Gene stable ID`)])==0)


ortho2 = merge(cynoanno2, ortho, all.x = T, sort = F)

ortho2[, seurat_gene := ifelse(`Gene name`%in% rownames(seurat_ape), `Gene name`, ifelse(`Gene stable ID` %in% rownames(seurat_ape), `Gene stable ID`, NA))]

genes_in_seurat = data.table(seurat_gene = rownames(seurat_ape))
ortho3 = merge(genes_in_seurat,ortho2,  by ="seurat_gene",  all.x = T, sort = F) %>% unique()

# ## how many genes do we win using orthologues instead of name matching ----
ortho3[`Gene name` != `Human gene name`, .(seurat_gene, `Gene name`, `Human gene name`)]
ortho3[is.na(`Gene name`)==F & seurat_gene != `Human gene name`, .(seurat_gene, `Gene name`, `Human gene name`)]
ortho3[is.na(`Gene name`)==F & seurat_gene != `Human gene name` & `Human gene name` %in% rownames(seurat_hum), .(seurat_gene, `Gene name`, `Human gene name`)]

# # check overlapping genes not in orthotable -----
qlist8  = venn4(ortho3$seurat_gene, rownames(seurat_ape),  rownames(seurat_hum),ortho3$`Human gene name`)
str(qlist8)
qlist8$q4 %>% sort() 

# add those
ortho3[, human_samename_only10x := ifelse(seurat_gene %in% rownames(seurat_hum) & is.na(`Human gene name`)==T & seurat_gene %nin% `Human gene name`,T, F)] # not those, already
ortho3[, uniqueN(seurat_gene),human_samename_only10x]

ortho3[human_samename_only10x==T ,.(seurat_gene, `Gene name`, `Human gene name`)]

ortho3[human_samename_only10x==T, `Human gene name` := seurat_gene]

qlist8b  = venn4(ortho3$seurat_gene, rownames(seurat_ape),  rownames(seurat_hum),ortho3$`Human gene name`)

# ## suggest 1:1 match ----
# oberserved n:m relation
ortho3[allDuplicatedEntries(seurat_gene) , .(seurat_gene, `Gene name`, `Human gene name`)]
ortho3[is.na(`Human gene name`)==F][allDuplicatedEntries(`Human gene name`) , .(seurat_gene, `Gene name`, `Human gene name`)]

# ## match via cluster-wise max expression ----

seurat_ape@assays

ape_meanexp_pre = AverageExpression(
  seurat_ape,
  assays = 'RNA',
  features = NULL,
  return.seurat = FALSE,
  group.by = "louvain_res_7",
  add.ident = NULL,
  layer = "data",
  verbose = TRUE
)

ape_meanexp =  as.data.table(ape_meanexp_pre$RNA,keep.rownames = T)
ape_meanexp[, max_across_cluster := max(.SD), rn]

### add to alliance_ortho
ortho3[, max_expression_cyno := ape_meanexp[match_hk(seurat_gene, ape_meanexp$rn), max_across_cluster]]

seurat_hum@assays

human_meanexp_pre = AverageExpression(
  seurat_hum,
  assays = 'RNA',
  features = NULL,
  return.seurat = FALSE,
  group.by = "louvain_res_7",
  add.ident = NULL,
  layer = "data",
  verbose = TRUE
)

human_meanexp =  as.data.table(human_meanexp_pre$RNA,keep.rownames = T)
human_meanexp[, max_across_cluster := max(.SD), rn]

### add to alliance_ortho
ortho3[, max_expression_human := human_meanexp[match_hk(`Human gene name`, human_meanexp$rn), max_across_cluster]]



qlist5 = venn2(rownames(seurat_ape), ortho3$seurat_gene)
stopifnot(all(rownames(seurat_ape) %in% ortho3$seurat_gene))


qlist5 = venn2(rownames(seurat_hum), ortho3$`Human gene name`)
qlist5$q3 %>% sort()

## assign best match to have 1:1
setorder(ortho3, seurat_gene, -max_expression_human,`Human gene name`, na.last = T)
ortho3[, ortho_1to1_cyno2human := duplicated(seurat_gene)==FALSE]
ortho3[allDuplicatedEntries(seurat_gene) , .(seurat_gene, `Gene name`, `Human gene name`, max_expression_cyno, max_expression_human, ortho_1to1_cyno2human)]


setorder(ortho3, `Human gene name`, -max_expression_cyno ,seurat_gene, na.last = T)
ortho3[is.na(`Human gene name`)==F, ortho_1to1_human2cyno := duplicated(`Human gene name`)==FALSE]
ortho3[is.na(`Human gene name`)==F][allDuplicatedEntries(`Human gene name`) , .(seurat_gene, `Gene name`, `Human gene name`, max_expression_cyno, max_expression_human, ortho_1to1_human2cyno)]


ortho3[,.N, .(ortho_1to1_human2cyno, ortho_1to1_cyno2human)]
ortho3_1to1 = ortho3[ortho_1to1_human2cyno==T & ortho_1to1_cyno2human==T]




qlist11 = venn2(ortho3_1to1$`Human gene name`, rownames(seurat_hum))


renamer = ortho3_1to1[`Human gene name` %in% qlist11$q1, .(seurat_gene, `Human gene name`)]
renamer[allDuplicatedEntries(seurat_gene)]
renamer[allDuplicatedEntries(`Human gene name`)]


# # save ----
fwrite(ortho3, here("analysisR/results/s0120_cyno2human_n_to_m_orthologues.csv.gz"))
fwrite(ortho3_1to1, here("analysisR/results/s0120_cyno2human_1to1_orthologues.csv.gz"))

fwrite(ortho3, here("analysisR/results_GIT//s0120_cyno2human_n_to_m_orthologues.csv.gz"))
fwrite(ortho3_1to1, here("analysisR/results_GIT//s0120_cyno2human_1to1_orthologues.csv.gz"))


fwrite(ortho3, here("analysisR/results/s0120_cyno2human_n_to_m_orthologues.csv.gz"))
fwrite(ortho3_1to1, here("analysisR/results/s0120_cyno2human_1to1_orthologues.csv.gz"))


saveRDS(seurat_ape, here("analysisR/results/s0120_seurat_QCed_cyno.rds"))
saveRDS(seurat_hum, here("analysisR/results/s0120_seurat_QCed_human.rds"))



# # finalize ----
finalizeSkript()