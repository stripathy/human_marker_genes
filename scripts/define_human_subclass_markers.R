library(dplyr)
library(Signac) # package used to find markers - no idea how good it is: 
library(parallel)

# read in pre-saved human snucseq gene expression file
human_exon_intron_cpm = readRDS(file = '~/allen_human/data/human_exon_intron_cpm.rds')
allen_human_data_dir = '~/allen_human/data-raw/allenHuman/'
allen_human_data_output_dir = '~/allen_human/data/'


humanGenes = read_csv(paste0(allen_human_data_dir, 'human_MTG_2018-06-14_genes-rows.csv'))
rownames(human_exon_intron_cpm) = humanGenes$gene

humanMetaJoined = readRDS( 
        file= paste0(allen_human_data_output_dir, 'humanMetaJoined.rds'))

# define human markers


allen_human_subclass_names = humanMetaJoined$new_subclass %>% levels
allen_human_subclass_names = allen_human_subclass_names[allen_human_subclass_names != 'NA.']

subclass_cluster_markers = mclapply(allen_human_subclass_names, function(curr_cluster_name){
  
  markers <- Signac::VeniceMarker(human_exon_intron_cpm[, !is.na(humanMetaJoined$new_subclass)], cluster = (humanMetaJoined[!is.na(humanMetaJoined$new_subclass), ]$new_subclass != curr_cluster_name), threshold = 10)
  new_cell_markers = markers %>% filter(Up.Down.score < -.99, Log2.fold.change < -5) %>% pull(Gene.Name) %>% as.character()
  return(new_cell_markers)
}, mc.cores = 20)
names(subclass_cluster_markers) = allen_human_subclass_names


saveRDS(subclass_cluster_markers, 
        file= paste0(allen_human_data_output_dir, 'derived_human_subclass_markers.rds'))


# test MGPs using aging dataset

# attempt MGP analysis
estimations =  mgpEstimate(exprData=mgp_gene_exp_mat_merged_trans,
                           genes=mouseMarkerGenes$Cortex,
                           geneColName='Gene.Symbol',
                           outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                           geneTransform =function(x){homologene::mouse2human(x)$humanGene}, # this is the default option for geneTransform
                           groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
                           seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                           removeMinority = TRUE)


# attempt MGP analysis
estimations_human_markers =  mgpEstimate(exprData=mgp_gene_exp_mat_merged_trans,
                           genes=subclass_cluster_markers,
                           geneColName='Gene.Symbol',
                           outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                           geneTransform =NULL, # this is the default option for geneTransform
                           groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
                           seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                           removeMinority = TRUE)

estimations$removedMarkerRatios


getMGPQualityMetrics = function(mgp_estimates){
  
  pc1_var_exp = lapply(names(mgp_estimates$trimmedPCAs), function(cell_type_name){
    
    curr_summary = mgp_estimates$trimmedPCAs[[cell_type_name]] %>% summary()
    pct_var_exp = curr_summary[["importance"]][2,1]
  }) %>% unlist
  names(pc1_var_exp) = names(mgp_estimates$trimmedPCAs)
  
  mgp_qual_df = data.frame(cell_type = names(mgp_estimates$trimmedPCAs), 
                           removed_marker_ratios = mgp_estimates$removedMarkerRatios)
  mgp_qual_df = cbind(mgp_qual_df, pc1_var_exp)
  return(mgp_qual_df)
}

mgp_qual_df %>% ggplot(aes(x = cell_type, y = removed_marker_ratios)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

mgp_qual_df %>% ggplot(aes(x = cell_type, y = pc1_var_exp)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

cell_type_list_ogan = c('Astrocyte', 'Endothelial', 'GabaPV', 'GabaRelnCalb', 'GabaVIPReln', 'Microglia',
                        'Microglia_activation', 'Microglia_deactivation', 'Oligo', 'OligoPrecursors', 'Pyramidal')
my_cell_type_list = c('Astrocyte', 'Endothelial', 'PVALB', 'SST', 'VIP', 'Microglia', NA, NA, 'Oligo', 'OPC', 'Pyramidal')


estimations_human_markers_quality = getMGPQualityMetrics(estimations_human_markers)
estimations_human_markers_quality$marker_source = 'human'

estimations_human_markers_quality$cell_type =  plyr::mapvalues(estimations_human_markers_quality$cell_type, 
                                                               my_cell_type_list, cell_type_list_ogan)

estimations_mouse_markers_quality = getMGPQualityMetrics(estimations)
estimations_mouse_markers_quality$marker_source = 'mouse'

combined_quality_df = rbind(estimations_mouse_markers_quality, estimations_human_markers_quality)
combined_quality_df$marker_source = factor(combined_quality_df$marker_source, levels = c('mouse', 'human'))

combined_quality_df %>% filter(!cell_type %in% c('Microglia_activation', 'Microglia_deactivation', 'LAMP5', 'PAX6'))


p1 = combined_quality_df %>% filter(!cell_type %in% c('Microglia_activation', 'Microglia_deactivation', 'LAMP5', 'PAX6')) %>% 
  ggplot(aes(x = cell_type, y = pc1_var_exp, fill = marker_source)) + 
  geom_bar(stat = "identity", position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = 'PC1 Variance Explained (ratio)')

save_plot(filename = '~/human_marker_genes/plots/mgp_qual_metrics_pc1.png', plot = p1, base_width = 8)

p2 = combined_quality_df %>% filter(!cell_type %in% c('Microglia_activation', 'Microglia_deactivation', 'LAMP5', 'PAX6')) %>% 
  ggplot(aes(x = cell_type, y = removed_marker_ratios, fill = marker_source)) + 
  geom_bar(stat = "identity", position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = 'Removed markers (ratio)')
save_plot(filename = '~/human_marker_genes/plots/mgp_qual_metrics_removed_markers.png', plot = p2, base_width = 8)



estimations$trimmedPCAs$GabaRelnCalb %>% summary()

plotMarkers(subclass_cluster_markers$PVALB, mouse_genes = F)

humanMetaNew = humanMeta
humanMetaNew$cluster_id = humanMetaNew$tree_order
humanMetaNew$cluster_label = humanMetaNew$cluster_name

lapp


