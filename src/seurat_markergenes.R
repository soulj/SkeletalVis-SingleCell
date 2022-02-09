library("optparse")

#Parse command line parameters
option_list <- list()

option_list$seuratObject <- make_option('--seuratObject', type='character' ,help="The path to the SeruatObject")
option_list$seuratOutput <- make_option('--seuratOutput', type='character', help="Seruat object with clusters annotated")
option_list$seuratDietOutput <- make_option('--seuratDietOutput', type='character', help="Slim seruat object with clusters annotated")
option_list$markerPlot <- make_option('--markerPlot', type='character', help="tSNE plot showing clusters labeled with marker genes")
option_list$diffExpTable <- make_option('--diffExpTable', type='character', help="tsv file with the differential expression between indentified clusters")


opt <- parse_args(OptionParser(option_list=option_list,description="Identify clusters and marker genes"))

#load needed libraries
library(ggplot2)
library(tidyverse)
library(Seurat)

#read in the seurat data
seu <- readRDS(opt$seuratObject)

#identify clusters of cells given a resolution
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu,resolution = 0.3)
seu <- FindVariableFeatures(seu, nfeatures = 3000)

#find the markers between each cluster and the others
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1)

#get the top markers for each cluster
markers.filt <- markers %>% group_by(cluster) %>% slice_min(n = 2, order_by = p_val_adj, with_ties = FALSE) %>% as.data.frame()

#format the top markers to plot
topMarkers <- markers.filt$gene
topMarkers <- split(topMarkers, ceiling(seq_along(topMarkers)/2))
topMarkers <- sapply(topMarkers,paste,collapse=",")

#write out the full list of markers for each cluster
write.table(markers, file=opt$diffExpTable, sep = "\t", row.names=F, quote = F, col.names = T)

#relabel the clusters
current.cluster.ids <- 0:(length(topMarkers)-1)
new.cluster.ids <- topMarkers
seu@active.ident <- plyr::mapvalues(x = seu@active.ident, from = current.cluster.ids, to = new.cluster.ids)

#visulise the clusters with the marker genes
g <- TSNEPlot(seu, label=TRUE, label.size=4)  + cowplot::theme_cowplot() + theme(legend.position = "none")
cowplot::save_plot(g, filename = opt$markerPlot, base_width=5, base_height=4, bg="white")


saveRDS(seu,opt$seuratOutput)

#create a low memory version for visualisation
seu <-DietSeurat(seu,counts = FALSE, dimreducs = c("tsne"), scale.data = FALSE)

saveRDS(seu,opt$seuratDietOutput)