convert_pathname_unix_windows = function(mypath, unixprefix = "/net/ifs1/san_projekte/projekte/genstat/", windowsprefix = "J:/", given = "unix") {
  # mypath = "/net/ifs1/san_projekte/projekte/genstat/02_projekte/2403_scRNAseq_cross_species_primate_human/establishing/scRNAseq_cross_species_primate_human_v01/../../scads/primates/reupload/data/run_nf/outdir_human_ensemble/human/cellranger/count//Human1_24hr_S5/outs/filtered_feature_bc_matrix.h5" 
  
  if(given == "unix") {
    mypath = str_replace_all(mypath, unixprefix, windowsprefix)
    
  } else if(given == "windows") {
    mypath = str_replace_all(mypath, windowsprefix, unixprefix)
  }
  mypath
}