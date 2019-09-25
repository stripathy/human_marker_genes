library(scrattch.vis)
library(tasic2016data)


options(stringsAsFactors = F)

anno <- tasic_2016_anno
anno <- anno[anno$primary_type_id > 0,]
data <- tasic_2016_rpkm[, anno$sample_name]

tasic_2016_genes = data %>% rownames()

data("mouseMarkerGenes")

use_genes = intersect(tasic_2016_genes, make.names(new_pvalb_markers))

# use_genes = new_pvalb_markers
data_df <- cbind(sample_name = colnames(data),
                 as.data.frame(t(data[use_genes,])))

# use_genes_plot = intersect(tasic_2016_genes, make.names(new_pvalb_markers))


group_dot_plot(data_df, 
                  anno, 
                  genes = use_genes, 
                  grouping = "primary_type", 
                  log_scale = FALSE,
                  font_size = 5,
                  rotate_counts = TRUE)

library(markerGeneProfile)

library(Signac)
library("Matrix")

library(tidyr)

class(pbmc.mat) 

anno = anno %>% separate(primary_type_label, into = c('subclass', 'marker_gene'), remove = F)

markers <- Signac::VeniceMarker(as(data, "dgCMatrix"), cluster = (anno$subclass != "L5b"), threshold = 10)


new_pvalb_markers = markers %>% filter(Up.Down.score < -.98, Log2.fold.change < -4) %>% pull(Gene.Name)
