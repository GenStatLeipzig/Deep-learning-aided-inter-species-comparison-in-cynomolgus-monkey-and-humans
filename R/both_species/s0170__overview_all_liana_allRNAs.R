rm(list = ls())
require(toolboxH)
require(here)

require(plotly)
require(ggthemes)
require(patchwork)
options(width = 322)

gzfiles_human_all  = dir(here("J:/02_projekte/2403_scRNAseq_cross_species_primate_human/github/scRNAseq_cross_species_primate_human_km/analysis/results/new_29_01_25/H26/"),pattern = "*iana_res.txt.gz", full.names = T)
gzfiles_human_all

gzfiles_cyno_all  = dir(here("J:/02_projekte/2403_scRNAseq_cross_species_primate_human/github/scRNAseq_cross_species_primate_human_km/analysis/results/new_29_01_25/M26/"),pattern = "*iana_res.txt.gz", full.names = T)
gzfiles_cyno_all


all = lapply(c(gzfiles_human_all, gzfiles_cyno_all) , function(myfile) {
  # myfile =  gzfiles[1]
  resi = fread(myfile)
  resi[, fn := myfile]
  resi

}) %>% rbindlist()


all[, csv := str_split(fn, "/") %>% sapply(last)]
all[, .N, csv]
all[, species := fifelse(grepl("^H",csv), "H.sapiens", "C.fascicularis")]
all[, .N,species]

all[, group1 := str_split(csv, "_", simplify = T)[,2] ]
all[,.N, group1]



hist(all$cellchat_pvals, breaks = 100)

all = moveColFront(all,c("species","group1"))
all[, log10magnitude_rank := -log10(magnitude_rank)]
all[, log10specificity_rank := -log10(specificity_rank)]
all[, log10specificity_rank := -log10(specificity_rank)]


require(ggplot2)
all[, n_targets := uniqueN(receptor_complex), .( group1, species, target)]
all[, n_sources := uniqueN(ligand_complex ), .( group1, species, source)]


all[,.N, .(specificity_rank<=0.05, magnitude_rank<=0.05)]

all2 = all[specificity_rank<=0.05 & magnitude_rank<=0.05]

 

all2_source = all2[ ,.(sum_observation = .N,
                       sum_magnitude = sum(log10magnitude_rank),
                       sum_specificity = sum(log10specificity_rank)), .( source, species, group1, n_cells_in_source_cluster, n_sources)]


all2_target = all2[ ,.(sum_observation = .N,
                       sum_magnitude = sum(log10magnitude_rank),
                       sum_specificity = sum(log10specificity_rank)
                       ), .( target, species, group1, n_cells_in_target_cluster, n_targets)]


# # plotting across celltypes ------------------------------------------------

## For which species should be ordered
species_ordered = "H.sapiens"

all2_source[,maxsum_magnitude_rel_source := max(sum_magnitude[species == species_ordered]/n_sources[species == species_ordered]) %>% as.numeric(), by = .(source)]
all2_target[,maxsum_magnitude_rel_target := max(sum_magnitude[species == species_ordered]/n_targets[species == species_ordered]) %>% as.numeric(), by = .(target)]

# # paperbilder ----

set.seed(16)
NogpaletteReihe <-  c("#CB769E", "#DE639A", "#A85C85", "#0081AF", "#4F6D7A", "#7C6A0A", "#368F8B", "#246A73", "#5CC1BC", "#62C370", "#F7C548", "#F97E44", "#FB3640", "#B7245C", "#0D3B66", "#3E2F5B", "#B2675E", "#644536")


NogpaletteReihe13 <-  c( "#0081AF", "#CB769E", "#DE639A",  "#A85C85","#368F8B",
                         "#5CC1BC","#B2675E","#B7245C", "#7C6A0A", "#644536",
                        "#62C370","#F7C548","#F97E44", "#FB3640",
                        "#4F6D7A", "#0D3B66", "#246A73", "#3E2F5B")


NogpaletteReihe14 <-  c( "#368F8B", alpha("#A85C85", 0.5), "#A85C85", "#DE639A",  "#368F8B",
                         "#5CC1BC","#B2675E","#B7245C", "#7C6A0A", "#644536",
                         "#62C370","#F97E44", "#FB3640",
                         "#4F6D7A", "#0D3B66", "#246A73", "#3E2F5B")

qlist2 = venn2(all2_target[species==species_ordered, target ],
               all2_target[species==species_ordered, target ], c(species_ordered, species_ordered))





p_targets = ggplot(all2_target[target %in% qlist2$q1], 
                     aes(reorder(target %>% str_replace_all("_", " ") %>% str_wrap(width = 17), -maxsum_magnitude_rel_target),
                         sum_magnitude/n_targets, 
                         fill = group1)) + 
  geom_col(position = position_dodge(width = 0.9, preserve = "single")) + 
  facet_wrap(~species, scales = "free_x") +
  scale_y_continuous(breaks = 0:100) +
  coord_flip() +
  theme_pander(base_size = 16) +
  xlab("") +
  ylab("Sum of normalized magnitude cell-cell interaction - \ntargets (Targets)") +
  scale_fill_manual(values = NogpaletteReihe14, 
                    guide = guide_legend(reverse = TRUE)) +
  ggtitle("Targets")


p_targets



p_sources = ggplot(all2_source[source %in% qlist2$q1], 
                  aes( reorder(source %>% str_replace_all("_", " ") %>% str_wrap(width = 17), -maxsum_magnitude_rel_source),
                       sum_magnitude/n_sources,
                       fill = group1)) +
  geom_col(position = position_dodge(width = 0.9, preserve = "single")) + 
  facet_wrap(~species, scales = "free_x")+
  scale_y_continuous(breaks = 0:100)+
  coord_flip()+
  theme_pander(base_size = 16) +
  xlab("")+
  ylab("Sum of normalized magnitude cell-cell interaction - \nsources (Sources)") +
  scale_fill_manual(values = NogpaletteReihe14, guide = guide_legend(reverse = TRUE) ) +
  ggtitle("Sources")
p_sources


p_count_all = p_targets / p_sources + plot_layout(guides = "collect") + plot_annotation(tag_suffix =  ")", tag_levels = "A", title = "All Sources-Receptor-Interactions", caption = "J:/02_projekte/2210_snRNA_Sympath_hi_maus_herz/scripts/s0170_magnitude_targets_Sources_grouped_allRNAs.pdf")  &
  theme(plot.tag = element_text(face = 'bold', size = 21), legend.position = "top") &
  labs(fill = "") &
  ylim(c(0,5))

p_count_all

pdf(here("analysisR/results/s0170_magnitude_targets_Sources_grouped_allRNAs.pdf"), width = 11, height = 11)
p_count_all
dev.off()
#
fwrite(all2_source, here("analysisR/results/s0170_magnitude_targets_grouped_all2_source_allRNAs.csv"))
fwrite(all2_target, here("analysisR/results/s0170_magnitude_targets_grouped_all2_target_allRNAs.csv"))

pdf(here("analysisR/results_GIT/s0170_magnitude_targets_Sources_grouped_allRNAs.pdf"), width = 11, height = 11)
p_count_all
dev.off()
#
fwrite(all2_source, here("analysisR/results_GIT/s0170_magnitude_targets_grouped_all2_source_allRNAs.csv"))
fwrite(all2_target, here("analysisR/results_GIT/s0170_magnitude_targets_grouped_all2_target_allRNAs.csv"))

#
finalizeSkript()
