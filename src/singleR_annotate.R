library("optparse")

#Parse command line parameters
option_list <- list()

option_list$seuratObject <- make_option('--seuratObject', type='character' ,help="The path to the seurat object rds file")
option_list$seuratOutput <- make_option('--seuratOutput', type='character', help="seurat object with cells automatically annotated")
option_list$cellTypeDatabase <- make_option('--cellTypeDatabase', type='character', help="reference dataset to use to annotate cells")
option_list$cellTypePlot <- make_option('--cellTypePlot', type='character' ,help="tSNE plot showing cell annotation")
option_list$cellTypeStats <- make_option('--cellTypeStats', type='character' ,help="numbers of each type of cell detected")
option_list$cellTypeMarkers <- make_option('--cellTypeMarkers', type='character' ,help="marker genes for each cell type")

opt <- parse_args(OptionParser(option_list=option_list,description="Annotate and visualise cell types"))

#load needed libraries
library(ggplot2)
library(tidyverse)
library(Seurat)
library(celldex)
library(SingleR)

#read in the Seurat object
seu <- readRDS(opt$seuratObject)

#get the reference database
if( opt$cellTypeDatabase =="HumanPrimaryCellAtlasData"){
  ref.se <- celldex::HumanPrimaryCellAtlasData()
} else {
  ref.se <- celldex::MouseRNAseqData()
}

#predict the cell type of each cell using the reference database
pred <- SingleR(test = as.SingleCellExperiment(seu), ref = ref.se, assay.type.test=1,
                     labels = ref.se$label.main)

#get the number of each predicted cell type
predictionTable <- data.frame(table(pred$labels))
colnames(predictionTable) <- c("CellType","Number")
write.table(predictionTable, file=opt$cellTypeStats, sep = "\t",row.names=F,quote = F,col.names = T)

#add the cell type predictions to the seurat metadata for later use
seu[["SingleR.labels"]] <- pred$labels

#visualise the data by predicted cell type label
g <- DimPlot(seu,group.by = "SingleR.labels")
cowplot::save_plot(g,filename = opt$cellTypePlot)


#change the identity of the cells to predicted cell type
Idents(seu) <- seu@meta.data$SingleR.labels

#get the marker genes for each predicted cell type
if(length(table(seu[["SingleR.labels"]]))>1){
  markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.table(markers, file=opt$cellTypeMarkers, sep = "\t",row.names=F,quote = F,col.names = T)
}

#save seurat object for later use
saveRDS(seu,opt$seuratOutput)