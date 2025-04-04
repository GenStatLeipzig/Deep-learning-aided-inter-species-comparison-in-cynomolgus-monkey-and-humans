#install.packages("data.table")
#install.packages("ggplot2")

require(data.table)
require(magrittr)
require(fdrtool)
require(ggplot2)

#human
stats = fread("/home/sc.uni-leipzig.de/mh89exih/compute/primates/MH84_c_R/H82_dge_tab_edgeR.csv")
#ggplot(stats, aes(PValue        )) + geom_histogram() + facet_grid(celltype ~comparison.vs.00hr, scales = "free")
etastats = stats[,fdrtool::fdrtool(PValue, statistic = "pvalue", plot = FALSE) %>% as.data.table(), .(celltype,comparison.vs.00hr)] %>% .[,c('celltype','comparison.vs.00hr', "statistic","param.cutoff","param.N.cens","param.eta0","param.eta0.SE" ),with = FALSE] %>% unique()
fwrite(etastats, "/home/sc.uni-leipzig.de/mh89exih/compute/primates/MH84_c_R/MH84_c_R_etahuman.csv")

#cyno
stats = fread("/home/sc.uni-leipzig.de/mh89exih/compute/primates/MH84_c_R/M82_dge_tab_edgeR.csv")
#ggplot(stats, aes(PValue        )) + geom_histogram() + facet_grid(celltype ~comparison.vs.00hr, scales = "free")
etastats = stats[,fdrtool::fdrtool(PValue, statistic = "pvalue", plot = FALSE) %>% as.data.table(), .(celltype,comparison.vs.00hr)] %>% .[,c('celltype','comparison.vs.00hr', "statistic","param.cutoff","param.N.cens","param.eta0","param.eta0.SE" ),with = FALSE] %>% unique()
fwrite(etastats, "/home/sc.uni-leipzig.de/mh89exih/compute/primates/MH84_c_R/MH84_c_R_etacyno.csv")

writeLines(capture.output(sessionInfo()), "/home/sc.uni-leipzig.de/mh89exih/compute/primates/MH84_c_R/MH84_c_R_session_info.txt")
#etastats
