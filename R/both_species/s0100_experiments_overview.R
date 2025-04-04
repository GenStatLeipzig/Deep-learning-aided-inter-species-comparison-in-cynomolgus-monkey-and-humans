rm(list = ls())

require(toolboxH)
require(here)
require(ggplot2)
require(ggthemes)
require(scales)
require(patchwork)
# File overview ----


# We have many files to import ----
human_allfiles_pre = dir(here("../../scads/primates/reupload/data/run_nf/outdir_human_ensemble/human/cellranger/count/"), full.names = T, recursive = T)
cyno_allfiles_pre = dir(here("../../scads/primates/reupload/data/run_nf/output_dir/cyno/cellranger/count/"), full.names = T, recursive = T)

allfiles =  c(human_allfiles_pre[grepl("filtered_feature_bc_matrix\\.h5",human_allfiles_pre)==T],
              cyno_allfiles_pre[grepl("filtered_feature_bc_matrix\\.h5",cyno_allfiles_pre)==T])

todo = data.table(filtered_h5_fn = allfiles)
todo[, experiment := str_split(filtered_h5_fn, "/") %>% sapply(., function(x) grep("_S", x, value = T)) ]
todo$experiment
todo[, tissue := "PBMC"]
todo[, .N,tissue]
todo[, species := ifelse(grepl("Human",experiment), "human", # Streptococcus pneumoniae 
                         ifelse(grepl("Cyno",experiment), "cyno" , experiment))]

todo[,.N, species]


todo[, individual := str_split(experiment, "_") %>% sapply(., "[", 1) %>% tolower()]
todo[,.N, individual]

todo[, timepoint_ori := str_split(experiment, "_") %>% sapply(., "[", 2)]
todo[,.N, timepoint_ori]

todo[,timepoint := ifelse(timepoint_ori =="TimeZero", "00hr", 
                           ifelse(timepoint_ori =="6hr", '06hr', timepoint_ori))]
todo[,.N, .(timepoint_ori, timepoint)] 


todo[, pairing := str_split(experiment, "_") %>% sapply(., "[", 3) %>% tolower]
todo[,.N, pairing]




# add raw files paths-----
rawfiles_dt =  data.table(raw_h5_fn  = c(human_allfiles_pre[grepl("raw_feature_bc_matrix\\.h5",human_allfiles_pre)==T],
                                         cyno_allfiles_pre[grepl("raw_feature_bc_matrix\\.h5",cyno_allfiles_pre)==T]))


rawfiles_dt[, experiment := str_split(raw_h5_fn, "/") %>% sapply(., function(x) grep("_S", x, value = T)) ]
rawfiles_dt$experiment
todo[, raw_h5_fn := rawfiles_dt[match_hk(todo$experiment, rawfiles_dt$experiment), raw_h5_fn]]


# ad numbering used in cell identifyer as suffix
setorder(todo,experiment)

# todo = cbind(suffix = 1:nrow(todo), todo)
todo[, cell_prefix := paste0(individual, "_", timepoint, "__")]
todo[, old_corefn :=  experiment]
todo = moveColFront(todo, c(  "experiment", "tissue", "species", 
                            "individual", "timepoint_ori", "timepoint",
                            "cell_prefix","old_corefn", "pairing", "filtered_h5_fn","raw_h5_fn"))


# experiment coherent with new names
todo[, experiment := paste0(individual,"_", timepoint,"_",pairing )]

# writexl::write_xlsx(list(todo = todo),here("analysisR/results/s0100_scRNA_experimental.xlsx") )
write.delim(todo,here("analysisR/results/s0100_scRNA_experimental.txt")) 

finalizeSkript()

