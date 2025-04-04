#' ---
#' title: "s0105_import_1_seurat.R"
#' output:
#'   html_document:
#'     keep_md: true
#'     toc: true
#'     number_sections: false
#'     toc_depth: 3
#'     toc_float: true
#'     code_folding: show
#' ---


# ## INITIALIZE ----

.libPaths()

require(Seurat)
require(toolboxH)
require(here)
require(ggplot2)
require(ggthemes)
require(scales)
require(patchwork)
require(stringr)

source(here("analysisR/scripts/R_helperfunctions.R"))
convert_pathname_unix_windows

# ## Experiment  overview ----


todo_pre = fread(here("analysisR/results/s0100_scRNA_experimental.txt"))


todo = todo_pre[is.na(raw_h5_fn)==F]

# check whether system is windows
if (Sys.info()[['sysname']] == "Windows") {
  todo$filtered_h5_fn = todo$filtered_h5_fn %>% convert_pathname_unix_windows()
  todo$raw_h5_fn = todo$raw_h5_fn %>% convert_pathname_unix_windows()
}

stopifnot(all(file.exists(todo$filtered_h5_fn)))
                          
stopifnot(all(file.exists(todo[is.na(raw_h5_fn)==F, raw_h5_fn ])))


todo_ids = todo$experiment %>% unique() %>% sort()#%>% .[1:2]
todo_ids


# # define mt-genes ----
cynoanno = fread(here("preprocessing/ensembl/mart_export.txt.gz"))
cynoanno[,.N, `Chromosome/scaffold name`]
mt_genes_cyno_pre = cynoanno[grepl("MT", `Chromosome/scaffold name`), ifelse(`Gene name`== "", `Gene stable ID`, `Gene name`)]
mt_genes_cyno = grep("ENSMFAG00", mt_genes_cyno_pre, invert = T, value = T)
mt_genes_cyno_tag = paste(paste0("^", mt_genes_cyno, "$"), collapse = "|")

mt_genes_cyno_tag

## orthologues
cynoortho = fread(here("preprocessing/ensembl/mart_export_orthologues.txt.gz"))
cynoortho[, chr_cyno := cynoanno[match_hk(cynoortho$`Gene stable ID`, cynoanno$`Gene stable ID`, makeunique = T, importcol = cynoanno$`Chromosome/scaffold name`), `Chromosome/scaffold name`]]
cynoortho[, name_cyno := cynoanno[match_hk(cynoortho$`Gene stable ID`, cynoanno$`Gene stable ID`, makeunique = T, importcol = cynoanno$`Gene name`), `Gene name`]]
cynoortho[, human_description := cynoanno[match_hk(cynoortho$`Gene stable ID`, cynoanno$`Gene stable ID`, makeunique = T, importcol = cynoanno$`Gene description`), `Gene description`]]

cynoortho[, samename := name_cyno == `Human gene name`]
notsamename_cyno   = cynoortho[samename == F & `Human gene name`!="" & name_cyno!=""]

notsamename_cyno[,.N, .(str_split(human_description, "\\[") %>% sapply(., "[", 1))][order(N)]


cynoortho_mt = cynoortho[chr_cyno =="MT" & name_cyno != ""]
cynoortho_mt

# # Import ----
all_obj = lapply(todo_ids, function(myid, mt_tag ="^MT-", rename_mt_table = cynoortho_mt) {
  # myid = todo$experiment[7]
  # mt_tag ="^MT-"
  # rename_mt_table = cynoortho_mt
  
  myrow = todo[experiment == myid]
    myfiltered_h5_fn = unique(myrow$filtered_h5_fn)
  stopifnot(length(myfiltered_h5_fn)==1)
  myobject = Read10X_h5(myfiltered_h5_fn, use.names = T)
  
  # look at the structure
  message("Found dimensions: ", paste(dim(myobject), collapse = " / "))
  
  old_dimnames = myobject@Dimnames[1] %>% unlist()
  if(any(old_dimnames %in% rename_mt_table$name_cyno)) {

  message("Renaming following MT genes from :",old_dimnames[match_hk(rename_mt_table$name_cyno,old_dimnames )]  %>% paste(., collapse = ", "))
  
  rownames(myobject)[match_hk(rename_mt_table$name_cyno,old_dimnames )] = rename_mt_table[match_hk(rownames(myobject)[match_hk(rename_mt_table$name_cyno,old_dimnames )], rename_mt_table$name_cyno), `Human gene name`]
   
  new_dimnames = myobject@Dimnames[1] %>% unlist()
  message("Renamed following MT genes to :",new_dimnames[match_hk(rename_mt_table$`Human gene name`,new_dimnames )]  %>% paste(., collapse = ", "))
}
  hh(myobject)
  
  # rename cell ids to allow unique names
  message("Renaming cell ids with prefix: ", myrow$cell_prefix)
  print(hh(myobject))
  colnames(myobject) = paste0(myrow$cell_prefix, colnames(myobject))
print(hh(myobject))

  # Setup Seurat object - filter for non-empty cells
  
  seurat_object_use= CreateSeuratObject(myobject, min.cells = 1, min.features = 1) #  more liberal than standard as i will filter min.features, i.e. cells expressen a certain a certain number of features later
  seurat_object_use$orig.ident %>% table()
  dim_after = dim(seurat_object_use) %>% huebsch()
  
  message("Importing seurtat for  ", myid,  " with dimensions " , paste0(dim_after, collapse = " x "), " genes x cells.")  
  seurat_object_use
  
 
    message("Using followin Mt-genes:\n", paste( grep(mt_tag, rownames(seurat_object_use), value=T, ignore.case = T), collapse = ", "))
  seurat_object_use[["pct.mito"]] <- PercentageFeatureSet(seurat_object_use, pattern = mt_tag) # The calculation here is uses simply the column sum of the matrix present in the counts slot
  p0_initQC =VlnPlot(seurat_object_use, features = c("nFeature_RNA", "nCount_RNA", "pct.mito"), ncol = 3, pt.size = 0) + plot_annotation(title = myrow$experiment, subtitle = myrow$group, caption = myfiltered_h5_fn) &  theme(axis.title.x = element_blank(),axis.text.x = element_blank()) 
  print(p0_initQC)
  
  
  
  vars2add = setdiff(names(todo), names(seurat_object_use@meta.data))
  vars2add
  # DefaultAssay(seurat_object_use)<-"RNA"
  
  for(i in vars2add){
    # i = vars2add[1]
    seurat_object_use[[i]] = myrow[[i]]  }
  
  seurat_object_use
  
}) #lapply


names(all_obj) = todo_ids


all_obj_merged = merge(all_obj[[1]],y = all_obj[2:length(all_obj)], merge.data = TRUE)
all_obj_merged
all_obj_merged@assays$RNA$counts.1 %>% hh()

anno = all_obj_merged@meta.data %>% as.data.table(keep.rownames = T)
anno[, prefix_found := str_split(rn, "__") %>% sapply("[", 1)]
anno[, .N, .(prefix_found, cell_prefix)]
anno[, stopifnot(identical(paste0(prefix_found, "__"), cell_prefix))]


# # look at mt genes per gene and individual----
all_obj_merged_joined = JoinLayers(all_obj_merged )
all_obj_merged_joined_mt = all_obj_merged_joined[grep("^MT-", rownames(all_obj_merged_joined)),]
mt_matrix = all_obj_merged_joined_mt@assays$RNA$counts
rownames(mt_matrix)
hh(mt_matrix) %>% as.matrix() %>% base::t()

mt_dt  = as.data.table(mt_matrix%>% as.matrix() %>% t(), keep.rownames = T)
mt_dt[, individual := str_split(rn, "_") %>% sapply("[", 1)]
mt_dtm = melt(mt_dt[, -"rn"], id.vars = "individual", variable.name = "gene", value.name = "counts")
mt_dtm[,sum(counts), .(gene, individual)] %>% dcast.data.table(gene ~ individual, value.var = "V1") %>% .[order(as.character(gene)), ] 
qlist2 = venn2(cynoortho_mt$`Human gene name`, rownames(all_obj_merged)) 
qlist2$q1
qlist2$q2


VlnPlot(all_obj_merged_joined, cynoortho_mt$`Human gene name`, split.by = 'species', ncol = 3)
jpeg(here("analysisR/results/s0105_vlnplot_mt-genes.jpg"), 7,14, units = "in", res = 300, quality = 100)
VlnPlot(all_obj_merged_joined, cynoortho_mt$`Human gene name`, split.by = 'species', ncol = 3)
dev.off()

anno[, .(counts = sum(nCount_RNA)), . (species, experiment)] %>% ggplot(., aes(x = experiment, y = counts, fill = species)) + geom_bar(stat = "identity") + 
   theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.5)) + ggtitle("Sum of counts per experiment") + xlab("Experiment") + ylab("Sum counts") + scale_y_continuous(labels = scales::comma) + theme(legend.position = "top")

# # SAVE----
saveRDS(all_obj_merged,file =  here("analysisR/results/s0105_all_obj_snRNA.RDS"))




# # finalize ----
finalizeSkript()


# .libPaths("~/rpackages/angmar/")
# setwd("/net/ifs1/san_projekte/projekte/genstat/02_projekte/2208_sympath_scRNA_LungeAtheroPneumoMaus/")
# rmarkdown::render(here::here("scripte/s201_4_create_seurat_dehasing.R"),encoding="UTF-8")

