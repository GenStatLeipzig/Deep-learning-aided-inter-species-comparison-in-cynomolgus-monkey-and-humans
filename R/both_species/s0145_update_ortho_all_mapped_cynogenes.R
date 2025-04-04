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

library(Seurat)
library(SummarizedExperiment)
library(SeuratObject)
library(toolboxH)
library(patchwork)



# # Read the AnnData as SingleCellExperiment object ----



# # Make seurat and visualise ----
# Check existing assays and keep X only, the most recent in python

ape_all = fread(here("data/MH108_scanpy_genes.csv"), header = T)
human_all = fread(here("data/MH108human_scanpy_genes.csv"), header = T)
# # Make a table including all seurat cyno genes annotated with orthologues ----
ortho = fread(here("preprocessing/ensembl/mart_export_orthologues.txt.gz"), na.strings = c("NA", "NULL", ""))
qlist1 = venn2(ortho$`Human gene name`, human_all$scanpy_gene)

cynoanno= fread(here("preprocessing/ensembl/mart_export.txt.gz"), na.strings = c("NA", "NULL", ""))
qlist2 = venn3(cynoanno$`Gene name`,cynoanno$`Gene stable ID`, ape_all$scanpy_gene)
str(qlist2)
sort(qlist2$q7, decreasing = T)

qlist2b = venn3(cynoanno$`Gene name` %>% toupper(),cynoanno$`Gene stable ID`, ape_all$scanpy_gene %>% toupper())
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

ortho2[, seurat_gene := ifelse(`Gene name`%in% ape_all$scanpy_gene, `Gene name`, ifelse(`Gene stable ID` %in% ape_all$scanpy_gene, `Gene stable ID`, NA))]

genes_in_seurat = data.table(seurat_gene = ape_all$scanpy_gene)
ortho3 = merge(genes_in_seurat,ortho2,  by ="seurat_gene",  all.x = T, sort = F) %>% unique()

# ## how many genes do we win using orthologues instead of name matching ----
ortho3[`Gene name` != `Human gene name`, .(seurat_gene, `Gene name`, `Human gene name`)]
ortho3[is.na(`Gene name`)==F & seurat_gene != `Human gene name`, .(seurat_gene, `Gene name`, `Human gene name`)]
 ortho3[is.na(`Gene name`)==F & seurat_gene != `Human gene name` & `Human gene name` %in% human_all$scanpy_gene, .(seurat_gene, `Gene name`, `Human gene name`)]

# # check overlapping genes not in orthotable -----
# qlist8  = venn4(ortho3$seurat_gene, ape_all$scanpy_gene,  human_all$scanpy_gene,ortho3$`Human gene name`)
# str(qlist8)
# qlist8$q4     %>% sort() 

# add those
ortho3[, human_samename_only10x := ifelse(seurat_gene %in% human_all$scanpy_gene & is.na(`Human gene name`)==T & seurat_gene %nin% `Human gene name`,T, F)] # not those, already
ortho3[, uniqueN(seurat_gene),human_samename_only10x]

ortho3[human_samename_only10x==T ,.(seurat_gene, `Gene name`, `Human gene name`)]

ortho3[human_samename_only10x==T, `Human gene name` := seurat_gene]

qlist8b  = venn4(ortho3$seurat_gene, ape_all$scanpy_gene,  human_all$scanpy_gene,ortho3$`Human gene name`)

qlist8c  = venn3(ortho3$seurat_gene, ape_all$scanpy_gene, ortho3[is.na(`Human gene name`)==F,seurat_gene], mylabels = c("Ape genes\nin orthotable", "All ape genes", "Ape genes\nwith human orthologue"))
ortho3
# 

# # save ----
fwrite(ortho3, here("analysisR/results/s0145_cyno2human_n_to_m_orthologues_unfilteredApes.csv.gz"))


fwrite(ortho3, here("analysisR/results_GIT//s0145_cyno2human_n_to_m_orthologues_unfilteredApes.csv.gz"))







# # finalize ----
finalizeSkript()