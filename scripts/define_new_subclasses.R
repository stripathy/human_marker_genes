# remap cell type names into a new column called new_subclass - which is a mix of level2, level3 and level4 classes

humanMetaJoined = humanMeta
humanMetaJoined$cluster = plyr::mapvalues(humanMetaJoined$cluster, old_cluster_names, new_cluster_names)

human_cluster_meta = read_excel(paste0('~/human_marker_genes/data-raw/allenHuman/', 'allen_human_nucseq_cluster_meta.xlsx'))
humanMetaJoined = left_join(humanMetaJoined, human_cluster_meta, by = 'cluster')

humanMetaJoined$subclass = factor(humanMetaJoined$level3, levels = unique(human_cluster_meta$level3))

humanMetaJoined$subclass_id = as.integer(humanMetaJoined$subclass)
humanMetaJoined$subclass_label = humanMetaJoined$subclass
humanMetaJoined$subclass_color = humanMetaJoined$level3_color

humanMetaJoined$new_subclass = as.character(humanMetaJoined$subclass)
# humanMetaJoined$new_subclass[str_detect(humanMetaJoined$cluster, 'PAX6')] = 'PAX6'
# humanMetaJoined$new_subclass[str_detect(humanMetaJoined$cluster, 'LAMP5')] = 'LAMP5'
humanMetaJoined$new_subclass[str_detect(humanMetaJoined$new_subclass, "LAMP5/PAX6/Other")] = make.names("LAMP5/PAX6/Other")
humanMetaJoined$new_subclass[str_detect(humanMetaJoined$cluster, 'OPC')] = 'OPC'
humanMetaJoined$new_subclass[str_detect(humanMetaJoined$cluster, 'Oligo')] = 'Oligo'
humanMetaJoined$new_subclass[str_detect(humanMetaJoined$cluster, 'Exc')] = 'Pyramidal'
humanMetaJoined$new_subclass = factor(humanMetaJoined$new_subclass, levels = unique(humanMetaJoined$new_subclass) %>% make.names(unique = T))

saveRDS(humanMetaJoined, 
        file= paste0(allen_human_data_output_dir, 'humanMetaJoined.rds'))