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


