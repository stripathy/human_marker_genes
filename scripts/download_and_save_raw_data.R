# download some allen human MTG data

dir.create('data-raw/allenHuman',showWarnings = FALSE)
download.file('http://celltypes.brain-map.org/api/v2/well_known_file_download/694416044',destfile = 'data-raw/allenHuman/humanExp.zip',method = 'wget')
unzip('~/allen_human/data-raw/allenHuman/humanExp.zip',exdir = '~/allen_human/data-raw/allenHuman')


### load human MTG nuclei data
allen_human_data_dir = '~/allen_human/data-raw/allenHuman/'
allen_human_data_output_dir = '~/allen_human/data/'

humanIntrons = read_csv(paste0(allen_human_data_dir, 'human_MTG_2018-06-14_intron-matrix.csv'))
humanExons = read_csv(paste0(allen_human_data_dir, 'human_MTG_2018-06-14_exon-matrix.csv'))
humanGenes = read_csv(paste0(allen_human_data_dir, 'human_MTG_2018-06-14_genes-rows.csv'))
rownames(humanExons) = humanGenes$gene
rownames(humanIntrons) = humanGenes$gene
humanMeta = read_csv(paste0(allen_human_data_dir, 'human_MTG_2018-06-14_samples-columns.csv'))
assertthat::assert_that(all(humanExons$X1 == humanIntrons$X1) & all(humanExons$X1 == humanGenes$entrez_id))
assertthat::assert_that(all(colnames(humanExons) == colnames(humanIntrons)))
assertthat::assert_that(all(colnames(humanExons)[-1] ==humanMeta$sample_name))
# We used log2-transformed CPM of intronic plus exonic reads for both datasets.
humanSum = humanExons[,-1]+humanIntrons[,-1]

human_exon_intron_cpm = cpm(humanSum, log = F)
human_exon_intron_cpm = as(human_exon_intron_cpm, "dgCMatrix")

human_exon_intron_log_cpm = cpm(humanSum, log = T, prior.count = .25)
human_exon_intron_log_cpm = as(human_exon_intron_log_cpm, "dgCMatrix")

human_exon_cpm = cpm(humanExons[,-1], log = F)

human_exon_log_cpm = cpm(humanExons[,-1], log = T, prior.count = .25)


saveRDS(human_exon_intron_cpm, 
        file= paste0(allen_human_data_output_dir, 'human_exon_intron_cpm.rds'))
saveRDS(human_exon_intron_log_cpm, 
        file= paste0(allen_human_data_output_dir, 'human_exon_intron_log_cpm.rds'))
saveRDS(human_exon_cpm, 
        file= paste0(allen_human_data_output_dir, 'human_exon_cpm.rds'))
saveRDS(human_exon_log_cpm, 
        file= paste0(allen_human_data_output_dir, 'human_exon_log_cpm.rds'))
