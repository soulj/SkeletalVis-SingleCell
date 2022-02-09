library("optparse")

#Parse command line parameters
option_list <- list()

option_list$seuratObject <- make_option('--seuratObject', type='character' ,help="The path to the seurat object RDS file")
option_list$seuratOutput <- make_option('--seuratOutput', type='character', help="seurat object with cell cycle effects removed")
option_list$preNormPlot <- make_option('--preNormPlot', type='character', help="tSNE plot showing cell cycle effects")
option_list$postNormPlot <- make_option('--postNormPlot', type='character' ,help="tSNE plot after cell cycle effect normalisation")

opt <- parse_args(OptionParser(option_list=option_list,description="Visualise and remove cell cycle effects given a SeruatObject"))

#load needed libraries
library(Matrix)
library(ggplot2)
library(DropletUtils)
library(Matrix)
library(tidyverse)
library(Seurat)
library(ggpointdensity)
library(scales)

#read in the seyrat object
seu <- readRDS(opt$seuratObject)

#get the cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Normalise the data and visualise the cells by predicted cell cycle stage
seu <- NormalizeData(seu) %>% ScaleData()
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seu <- FindVariableFeatures(seu, nfeatures = 3000)
seu <- RunPCA(seu, verbose = FALSE, npcs = 20)
seu <-  RunTSNE(seu, dims = 1:10,features = c(s.genes, g2m.genes),check_duplicates = FALSE)
g <- DimPlot(seu)
cowplot::save_plot(g,filename = opt$preNormPlot)

#remove unwanted variation due to cell cycle or mitochondrial content
seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score","percent.mt"))

#visualise the new data
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu,resolution = 0.5)
seu <- FindVariableFeatures(seu, nfeatures = 3000)
seu <- RunPCA(seu, verbose = FALSE, npcs = 20)
seu <- RunTSNE(seu, dims = 1:10,features = c(s.genes, g2m.genes),check_duplicates = FALSE)


g <- TSNEPlot(seu,group.by="Phase",cols=hue_pal()(3)[c(2,3,1)])
cowplot::save_plot(g,filename = opt$postNormPlot)

saveRDS(seu,opt$seuratOutput)