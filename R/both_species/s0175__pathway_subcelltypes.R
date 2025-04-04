require(toolboxH)
require(here)
require(ggplot2)
require(plotly)
require(Seurat)
require(ggthemes)

# Load
seurat  = readRDS(here("analysisR/results/s0140_seurat2.rds"))



# Viz
seurat$species_timepoint = paste(seurat$species, seurat$timepoint, sep = "_")

Idents(seurat) %>% table()
DimPlot(seurat[,seurat$cluster_azimut1_5_scanvi_nkuni == "CD8 T"], label = TRUE, split.by = "species_timepoint",reduction =     "umap.scvi", ncol = 3)


# Diff exp cell-wise
seurat2 = JoinLayers(seurat)

seurat2$species %>% table()

cluster10.markers_cyno <- FindMarkers(seurat2[,seurat2$cluster_azimut1_5_scanvi_nkuni == "CD8 T" & seurat2$species == "cyno"], ident.1 = 10)  %>% data.table(keep.rownames = T)
cluster10.markers_cyno[, pvals_adj_0.05 := p_val_adj <=0.05]
min_expressed = 0.2
cluster10.markers_cyno[, expressed := fifelse((avg_log2FC >0 & pct.1 > min_expressed) | 
                                                (avg_log2FC <0 & pct.2 > min_expressed), "good", "low")]

cluster10.markers_human<- FindMarkers(seurat2[,seurat2$cluster_azimut1_5_scanvi_nkuni == "CD8 T" & seurat2$species == "human"], ident.1 = 10)  %>% data.table(keep.rownames = T)
cluster10.markers_human[, pvals_adj_0.05 := p_val_adj <=0.05]
cluster10.markers_human[, expressed := fifelse((avg_log2FC >0 & pct.1 > min_expressed) | 
                                                (avg_log2FC <0 & pct.2 > min_expressed), "good", "low")]

cluster10.markers = merge.data.table(cluster10.markers_human, cluster10.markers_cyno, by = c("rn"), suffixes = c("_human", "_cyno"), all = T)
cluster10.markers

p2 = ggplot(cluster10.markers, aes(x = avg_log2FC_human     , y = avg_log2FC_cyno ,label = rn, col= paste(pvals_adj_0.05_cyno , pvals_adj_0.05_human))) + geom_point() + geom_smooth(method = "lm") + geom_abline(lty = 3) + theme_minimal()
p2
ggplotly(p2)

p2 %+% cluster10.markers[expressed_cyno == "good" & expressed_human =="good"]


p2b = ggplot(cluster10.markers, aes(x = sign(avg_log2FC_human)*-log10(p_val_human)     , y = sign(avg_log2FC_cyno)*-log10(p_val_cyno) ,label = rn, col= paste(pvals_adj_0.05_cyno , pvals_adj_0.05_human))) + geom_point() + geom_smooth(method = "lm") + geom_abline(lty = 3) + theme_minimal()
p2b

p2b %+% cluster10.markers[expressed_cyno == "good" & expressed_human =="good"]







## pw enrich
dt = cluster10.markers[expressed_cyno == "good" & expressed_human =="good"]
library(clusterProfiler)
library(ReactomePA)
library(ggridges)
library(data.table)
library(org.Hs.eg.db)

# Prepare ranked gene lists
dt[, entrez := bitr(dt$rn, fromType="SYMBOL", toType="ENTREZID", 
                    OrgDb="org.Hs.eg.db", drop=F)$ENTREZ]

dt[, SYMBOLfromEntrez := bitr(dt$rn, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = F)$SYMBOL]
dt[,stopifnot(identical(rn, SYMBOLfromEntrez))]
dt$SYMBOLfromEntrez = NULL

# Create and sort ranked lists for both species
human_ranks <- dt[!is.na(entrez), setNames(avg_log2FC_human, entrez)] %>% sort(decreasing=TRUE)
cyno_ranks <- dt[!is.na(entrez), setNames(avg_log2FC_cyno, entrez)] %>% sort(decreasing=TRUE)

# Run GSEA for both species


# i use log2FC for nice interpretability
dt[, plot(avg_log2FC_human, -log10(p_val_human), col = ifelse(pvals_adj_0.05_human, "red", "black"), pch = 19, main = "Human", xlab = "log2FC", ylab = "-log10(p-value)")]

dt[, plot(avg_log2FC_cyno, -log10(p_val_cyno), col = ifelse(pvals_adj_0.05_cyno, "red", "black"), pch = 19, main = "CYNO", xlab = "log2FC", ylab = "-log10(p-value)")]


human_gsea <- gsePathway(human_ranks, pvalueCutoff=1, pAdjustMethod="BH", eps=0, by="fgsea")
cyno_gsea <- gsePathway(cyno_ranks, pvalueCutoff=1, pAdjustMethod="BH", eps=0, by="fgsea")

# Get top 10 pathways for each species
maxpath <- 10

human_gsea_dt[, first3 :=strsplit(core_enrichment, "/")[[1]][1:1] %>%   paste(collapse="/") , core_enrichment]
human_gsea_dt2 = human_gsea_dt[duplicated(first3)==F][order(-abs(NES))][1:maxpath]


human_top <- setReadable(human_gsea, OrgDb=org.Hs.eg.db, keyType="ENTREZID") %>% 
  as.data.table() %>% 
  .[, firstgenes :=strsplit(core_enrichment, "/")[[1]][1:1] %>%   paste(collapse="/") , core_enrichment] %>% .[duplicated(firstgenes)==F] %>% 
  .[order(-abs(NES))] %>% 
  .[1:maxpath, ]

cyno_top <- setReadable(cyno_gsea, OrgDb=org.Hs.eg.db, keyType="ENTREZID") %>% 
  as.data.table() %>% 
  .[, firstgenes :=strsplit(core_enrichment, "/")[[1]][1:1] %>%   paste(collapse="/") , core_enrichment] %>% .[duplicated(firstgenes)==F] %>% 
  .[order(-abs(NES))] %>% 
  .[1:maxpath, ]

# Combine unique pathways
combined_pathways <- unique(c(human_top$Description, cyno_top$Description))

# Filter GSEA results for these pathways
human_gsea_filt <- human_gsea
human_gsea_filt@result <- human_gsea_filt@result[human_gsea_filt$Description %in% combined_pathways,]
cyno_gsea_filt <- cyno_gsea
cyno_gsea_filt@result <- cyno_gsea_filt@result[cyno_gsea_filt$Description %in% combined_pathways,]

# Create combined ridge plot data
joint_data <- rbind(
  ridgeplot(human_gsea_filt)$data %>% data.table() %>% .[, species := "human"],
  ridgeplot(cyno_gsea_filt)$data %>% data.table() %>% .[, species := "cyno"]
)




## add NES
nesinfo = rbind(human_gsea_filt@result %>% data.table() %>% .[, species := "human"],
                cyno_gsea_filt@result %>% data.table() %>% .[, species := "cyno"]) %>% 
  .[Description %in% combined_pathways, .(Description, NES, species)] %>% 
  dcast.data.table(Description ~ species, value.var = "NES") %>% 
  .[order(sign(human), -abs(human), -abs(cyno))] %>%
  .[,.(Description,Description2= paste0(Description, " (H:", format(round(human, 1), nsmall=1), " | C:", format(round(cyno, 1), nsmall=1), ")") %>% str_wrap(55))]

nesinfo

# Create final plot
joint_data[, Description2 := nesinfo[match_hk(joint_data$category, nesinfo$Description), Description2]]

joint_data[, Description2 := factor(Description2, levels=nesinfo$Description2, ordered=TRUE)]
p_ridge_combined_bothcentric <- ggplot(joint_data, 
                                       aes(x=value, y=Description2, fill=species)) +
  geom_density_ridges(alpha=0.6) +
  theme_ridges(font_size = 16) +
  scale_fill_manual(values=c("blue", "red")) +
  theme(legend.position=c(0.2, 0.8)) +
  labs(title="Cluster 10 vs. others (Top 10 cyno & 10 human)") +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x="log2 Foldchange core pathway genes",
       y = "",
       fill = "") 


p_ridge_combined_bothcentric
pdf(here("analysisR/results_GIT/s0175_cluster10_cd8_vs_others_pathway.pdf"), 10.5,7)
p_ridge_combined_bothcentric
dev.off()
