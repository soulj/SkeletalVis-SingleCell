library("optparse")

#Parse command line parameters
option_list <- list()

option_list$seuratObject <- make_option('--seuratObjects', type='character', default = NULL,help="Integrated seurat object")
option_list$upRegulated <- make_option('--upRegulated', type='character', help="Differential expression from upregulated neighbourhoods")
option_list$downRegulated <- make_option('--downRegulated', type='character', help="Differential expression from downregulated neighbourhoods")
option_list$neighbourhoods <- make_option('--neighbourhoods', type='character', help="Plot of the differentially abundant neighbourhoods")

opt <- parse_args(OptionParser(option_list=option_list))

#load needed libraries
library(Matrix)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(cowplot)
library(miloR)
library(scater)
library(scran)

seu <- readRDS(opt$seuratObject)

sce <- as.SingleCellExperiment(seu)

sce_milo <- Milo(sce)

traj_milo <- buildGraph(sce_milo, k = 30, d = 30)
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 30, d=30, refined = TRUE)

plotNhoodSizeHist(traj_milo)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="orig.ident")

traj_design <- data.frame(colData(traj_milo))[,c("orig.ident", "Condition")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident
traj_milo <- calcNhoodDistance(traj_milo, d=30)
da_results <- testNhoods(traj_milo, design = ~ Condition, design.df = traj_design)

traj_milo <- buildNhoodGraph(traj_milo)

g <- plotUMAP(traj_milo) + plotNhoodGraphDA(traj_milo, da_results, alpha=0.1)
save_plot(g,file=opt$neighbourhoods,bg="white")

# ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
# 
# ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
#   geom_point() +
#   geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

# Add log normalized count to Milo object
da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC > 0)
da_nhood_markers <- findNhoodGroupMarkers(traj_milo, da_results,
                                          aggregate.samples = TRUE, sample_col = "orig.ident")
write.table(da_nhood_markers, file=opt$upRegulated, sep = "\t",row.names=F,quote = F,col.names = T)


da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
da_nhood_markers <- findNhoodGroupMarkers(traj_milo, da_results,
                                          aggregate.samples = TRUE, sample_col = "orig.ident")

write.table(da_nhood_markers, file=opt$downRegulated, sep = "\t",row.names=F,quote = F,col.names = T)

