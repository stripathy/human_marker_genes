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



# try using ogan's MGP stuff to find markers

human_exon_intron_log_cpm = readRDS('~/allen_human/data/human_exon_intron_log_cpm.rds')
allenGenes = read_csv('~/allen_human/data-raw/allenHuman/human_MTG_2018-06-14_genes-rows.csv')


library(markerGeneProfile)


humanMetaJoinedUse = humanMetaJoined %>% filter(new_subclass == 'Endothelial' | new_subclass == 'Microglia')

humanMetaJoinedUseFilt = humanMetaJoined %>% filter(!is.na(new_subclass),
                                                    !(new_subclass == 'Endothelial'), 
                                                    !(new_subclass == 'Microglia')
) %>%
  group_by(new_subclass) %>%
  sample_n(100, replace = F)

humanMetaJoinedUse = bind_rows(humanMetaJoinedUse, humanMetaJoinedUseFilt)


use_single_cell_ids = humanMetaJoinedUse %>% pull(sample_name) 
allen_sce = human_exon_intron_log_cpm[,use_single_cell_ids] %>% as.data.frame()

# need to load allen gene data
rownames(allen_sce) = allenGenes$gene

# allen_sce = human_exon_intron_log_cpm[,use_single_cell_ids] %>% as.matrix() %>% as.data.frame()
allen_sce %<>% tibble::rownames_to_column(var = 'Gene.Symbol')
# allen_sce %<>% mutate(Gene.Symbol = make.names(Gene.Symbol, unique = T))

# allen_sce %<>% filter(Gene.Symbol %in% c('SST', 'PVALB'))
# allen_sce = allen_sce %>% select(Gene.Symbol)

humanMetaMarkerSelect = humanMetaJoinedUse %>% select(sample_name, new_subclass) %>% mutate(new_subclass = as.character(new_subclass))
humanMetaMarkerSelect$replicate = 1:nrow(humanMetaMarkerSelect)

for (i in 1:10){
  markerCandidates(design = humanMetaMarkerSelect, # the design file
                   expression = allen_sce, # expression file 
                   outLoc = file.path('README_files/Rotation',i), # output directory
                   groupNames = 'new_subclass', # name of the column with cell types. can be a vector
                   regionNames = NULL, # name of the column with brain regions. leave NULL if no region seperation is desired
                   # PMID = 'PMID', # name of the column with study identifiers
                   sampleName = 'sample_name', # name of the column with sample names
                   replicates = 'replicate', # name of the column with replicates
                   foldChangeThresh = 10, # threshold of fold change for gene selection (default is 10)
                   minimumExpression = 3, # minimum expression level that a gene can be considered a marker gene (default is 8)
                   background = 0, # background level of expression (default is 6)
                   # regionHierarchy = mpg_sampleRegionHiearchy, # hierarchy of brain regions to be used
                   geneID = 'Gene.Symbol', # column name with with gene idenditifers
                   cores = 12, # number of cores to use in parallelization 
                   rotate = 0.33,
                   seed = i
  )
}

rotateSelect(rotationOut='README_files/Rotation',
             rotSelOut='README_files/RotSel',
             cores = 8,
             foldChange = 1 # this is a fixed fold change threshold that ignores some leniency that comes from markerCandidates. setting it to 1 makes it irrelevant
)

pickMarkers('README_files/RotSel/All_CellType/',rotationThresh = 0.95)




# human_sc_expr = 
markerCandidates(design = humanMetaMarkerSelect, # the design file
                 expression = allen_sce, # expression file 
                 outLoc = 'README_files/quickSelection', # output directory
                 groupNames = 'new_subclass', # name of the column with cell types. can be a vector
                 regionNames = NULL, # name of the column with brain regions. leave NULL if no region seperation is desired
                 sampleName = 'sample_name', # name of the column with sample names
                 replicates = 'replicate', # name of the column with replicates
                 foldChangeThresh = 10, # threshold of fold change for gene selection (default is 10)
                 minimumExpression = 1, # minimum expression level that a gene can be considered a marker gene (default is 8)
                 background = 0, # background level of expression (default is 6)
                 geneID = 'Gene.Symbol', # column name with with gene idenditifers
                 cores = 20 # number of cores to use in parallelization 
)


data(mgp_sampleProfiles)
data(mgp_sampleProfilesMeta)

markerCandidates(design = mgp_sampleProfilesMeta, # the design file
                 expression = mgp_sampleProfiles, # expression file 
                 outLoc = 'README_files/quickSelection', # output directory
                 groupNames = 'CellType', # name of the column with cell types. can be a vector
                 regionNames = NULL, # name of the column with brain regions. leave NULL if no region seperation is desired
                 PMID = 'PMID', # name of the column with study identifiers
                 sampleName = 'sampleName', # name of the column with sample names
                 replicates = 'replicate', # name of the column with replicates
                 foldChangeThresh = 10, # threshold of fold change for gene selection (default is 10)
                 minimumExpression = 8, # minimum expression level that a gene can be considered a marker gene (default is 8)
                 background = 6, # background level of expression (default is 6)
                 geneID = 'Gene.Symbol', # column name with with gene idenditifers
                 cores = 8 # number of cores to use in parallelization 
)

read.table('README_files/quickSelection/CellType/Cell C') %>% knitr::kable()

human_markers_quick_sel = pickMarkers('README_files/quickSelection/new_subclass/',
            foldChange = 1,  # this is a fixed fold change threshold that ignores some leniency that comes from markerCandidates. setting it to 1 makes it irrelevant
            silhouette = 0.5)

base_triplab_dir = '/external/rprshnas01/netdata_kcni/stlab/marker_genes/'
saveRDS(human_markers_quick_sel, 
        file= paste0(base_triplab_dir, 'human_markers_quick_sel.rds'))



