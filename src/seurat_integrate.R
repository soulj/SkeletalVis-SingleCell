library("optparse")

#Parse command line parameters
option_list <- list()

option_list$seuratObjects <- make_option('--seuratObjects', type='character', default = NULL,help="The path to seurat object RDS files to be integrated")
option_list$seuratOutput <- make_option('--seuratOutput', type='character', help="Output integrated Seruat object")
option_list$tSNEPlot <- make_option('--tSNEPlot', type='character', help="tSNEPlot coloured by sample id")

opt <- parse_args(OptionParser(option_list=option_list))

#load needed libraries
library(Matrix)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(cowplot)

#get the seurat objects to integrate
seu.list <- lapply(unlist(strsplit(opt$seuratObjects," ")),readRDS)

#normalise all the samples
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

#find the integration features
features <- SelectIntegrationFeatures(object.list = seu.list)
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features, reduction = "rpca",
                                  dims = 1:30)

#integrate the data
integrated <- IntegrateData(anchorset = anchors, dims = 1:30,k.weight = 20)

integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:30)

#visualise the integrated data
g <- DimPlot(integrated, group.by = "orig.ident")
cowplot::save_plot(g, filename = opt$tSNEPlot, base_width=5, base_height=4, bg="white")

saveRDS(integrated,opt$seuratOutput)
