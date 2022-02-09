library("optparse")

#Parse command line parameters
option_list <- list()

option_list$seuratObject <- make_option('--seuratObject', type='character' ,help="The path to the seurat object rds file")
option_list$groupBy <- make_option('--groupBy', type='character', help="Label to group the cells by")
option_list$interactions <- make_option('--interactions', type='character', help="Results table giving the ligand-receptor interactions for each celltype/cluster")

opt <- parse_args(OptionParser(option_list=option_list,description="Identify cell crosstalk with liana"))


#BiocManager::install("OmnipathR")
#devtools::install_github('saezlab/liana')
library(liana)
library(tidyverse)
library(Seurat)

seu <- readRDS(opt$seuratObject)

Idents(seu) <- opt$groupBy

liana_data <- liana_wrap(seu)
liana_data <- liana_data %>%
  liana_aggregate()

write.table(liana_data, file=opt$interactions, sep = "\t",row.names=FALSE,quote = FALSE,col.names = TRUE)
