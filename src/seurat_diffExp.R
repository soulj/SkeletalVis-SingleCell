library("optparse")

#Parse command line parameters
option_list <- list()

option_list$seuratObject <- make_option('--seuratObject', type='character' ,help="The path to the SeruatObject")
option_list$comparisonsTable <- make_option('--comparisonsTable', type='character', help="tsv file with the differential expression between indentified clusters")
option_list$diffExpTableClusters <- make_option('--diffExpTableClusters', type='character', help="tsv file with the differential expression within clusters between conditions")
option_list$diffExpTableCellTypes <- make_option('--diffExpTableCellTypes', type='character', help="tsv file with the differential expression in cell types between conditions")


opt <- parse_args(OptionParser(option_list=option_list,description="Identify clusters and marker genes"))

#load needed libraries
library(ggplot2)
library(tidyverse)
library(Seurat)

seu <- readRDS(opt$seuratObject)

seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu,resolution = 0.3)
seu <- FindVariableFeatures(seu, nfeatures = 3000)

#get the gene markers for each cluster between conditions
seu$celltype.condition <- paste(Idents(seu), seu$Condition, sep = "_")
seu$celltype <- Idents(seu)
Idents(seu) <- "celltype.stim"
markers <- FindAllMarkers(seu, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers, n = 15)

write.table(markers, file=opt$diffExpTableClusters, sep = "\t",row.names=F,quote = F,col.names = T)
