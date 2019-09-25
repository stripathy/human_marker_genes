library(readr)
library(edgeR)

# download reference files from Allen
dir.create('data-raw/allenMouse',showWarnings = T)
download.file('http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985',destfile = 'data-raw/allenMouse/mouseExp.zip')
unzip('data-raw/allenMouse/mouseExp.zip',exdir = 'data-raw/allenMouse')


### load mouse visual cortex dissociated cell data

# mouseIntrons = read_csv('data-raw/allenMouse/mouse_VISp_2018-06-14_intron-matrix.csv')
mouseExons = read_csv('data-raw/allenMouse/mouse_VISp_2018-06-14_exon-matrix.csv')
mouseGenes = read_csv('data-raw/allenMouse/mouse_VISp_2018-06-14_genes-rows.csv')
mouseMeta = read_csv('data-raw/allenMouse/mouse_VISp_2018-06-14_samples-columns.csv')
assertthat::assert_that(all(mouseExons$X1 == mouseIntrons$X1) & all(mouseExons$X1 == mouseGenes$entrez_id))
assertthat::assert_that(all(colnames(mouseExons) == colnames(mouseIntrons)))
assertthat::assert_that(all(colnames(mouseExons)[-1] == mouseMeta$sample_name))
# We used log2-transformed CPM of intronic plus exonic reads for both datasets.

# mouseSum = mouseExons[,-1] #+mouseIntrons[,-1]
mouseCPM = cpm(mouseExons[,-1])
# mouseCPMLog = log2(mouseCPM +1)
rownames(mouseCPM) = mouseGenes$gene_symbol

saveRDS(mouseCPM,file= 'data-raw/allenMouse/mouseAll.rds')


mouse_cpm_sparse = as(mouseCPM, "dgCMatrix")

# anno = anno %>% separate(primary_type_label, into = c('subclass', 'marker_gene'), remove = F)

markers <- Signac::VeniceMarker(mouse_cpm_sparse, cluster = (mouseMeta$cluster != "Pvalb Gabrg1"), threshold = 10)


new_cell_markers = markers %>% filter(Up.Down.score < -.98, Log2.fold.change < -5) %>% pull(Gene.Name)



use_genes = intersect(tasic_2016_genes, make.names(new_cell_markers))

# use_genes = new_cell_markers
data_df <- cbind(sample_name = colnames(data),
                 as.data.frame(t(data[use_genes,])))

# use_genes_plot = intersect(tasic_2016_genes, make.names(new_cell_markers))


group_dot_plot(data_df, 
               anno, 
               genes = use_genes, 
               grouping = "primary_type", 
               log_scale = FALSE,
               font_size = 5,
               rotate_counts = TRUE)