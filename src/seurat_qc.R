library("optparse")

#Parse command line parameters
option_list <- list()


option_list$tenx <- make_option('--tenx', type='character', default = NULL,help="The path to the 10x file from kallistobus")
option_list$cellranger <- make_option('--cellranger', type='character', default=NULL, help="The path to the cellranger directory")
option_list$t2g <- make_option('--t2g', type='character', help="The path to the transcript to gene mapping file from kallistobus")
option_list$sampleName <- make_option('--sampleName', type='character', default = NULL,help="The path to the 10x file from kallistobus")
option_list$metaData <- make_option('--metaData', type='character', default = NULL,help="MetaData file linking the sample name to experimental condition")
option_list$nFeatureMin <- make_option('--nFeatureMin', type='double', default =200 ,help="Minimum number of genes for cell filtering")
option_list$nFeatureMax <- make_option('--nFeatureMax', type='double', default =7000 ,help="Maximum number of genes for cell filtering")
option_list$percentMTMax <- make_option('--percentMTMax', type='double', default =10 ,help="Maximum cell mitochondrial content for filtering")
option_list$seuratOutput <- make_option('--seuratOutput', type='character', help="Seruat object with empty droplets removed")
option_list$QCPlot <- make_option('--QCPlot', type='character', help="QC plots showing droplets, mito,mRNA and gene stats")

opt <- parse_args(OptionParser(option_list=option_list))

#load needed libraries
library(Matrix)
library(ggplot2)
library(DropletUtils)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(cowplot)


#function to read in cell x gene matrix from kallisto-bustools and annotate with barcodes.
read_count_output <- function(dir, name) {

  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")

  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  closeAllConnections()
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

#TODO - test if both file options are null or not
metaDataTable <- read.delim(opt$metaData)

if(!is.null(opt$tenx)){

#read in the transcript to gene map
tr2g <- read_tsv(opt$t2g, col_names = c("transcript", "gene", "gene_name"))

#read in the cell x gene matrix
res_mat <- read_count_output(opt$tenx,"cells_x_genes")

#remove any empty rows
res_mat <- res_mat[Matrix::rowSums(res_mat) > 0,]

bcr <- barcodeRanks(res_mat)
knee <- metadata(bcr)$knee
inflection <- metadata(bcr)$inflection
bcr <- as.data.frame(bcr)

# Making a plot.

barcodePlot <- ggplot(bcr,aes(rank,total)) + geom_line()+ xlab("Rank") + ylab("Total") +
   scale_y_continuous(trans='log10') +  scale_x_continuous(trans='log10') +
  geom_hline(yintercept=knee,linetype="dotted") +
  geom_hline(yintercept=inflection,linetype="dotted") +
  annotate(
    "text", label = "Knee",
    x = 5, y = knee, size = 6, colour = "black"  ) +
  annotate(
    "text", label = "Inflection",
    x = 5, y = inflection, size = 6, colour = "black"  ) +
  cowplot::theme_cowplot(font_size = 20)


#reproducibly infer which droplets contain cells
set.seed(42)
e.out <- emptyDrops(res_mat)

#filter droplets
is.cell <- which(e.out$FDR <= 0.01)
res_mat <- res_mat[,is.cell]
res_mat <- res_mat[Matrix::rowSums(res_mat) > 0,]

e.out[ is.na(e.out$FDR),"FDR"] <- 1
is.cell <- e.out$FDR <= 0.01

cellPlot  <- ggplot(as.data.frame(e.out),aes(Total,-LogProb,color=is.cell)) + geom_point()+
  xlab("Total UMI count") + ylab ("-Log Probability") +
  scale_y_continuous(trans='log10') +  scale_x_continuous(trans='log10') +
  cowplot::theme_cowplot(font_size = 20)
  
# Convert from Ensembl gene ID to gene symbol
rownames(res_mat) <- tr2g$gene_name[match(rownames(res_mat), tr2g$gene)]

} else if(!is.null(opt$cellranger)){

res_mat <- Read10X(opt$cellranger)

}

#create an initial Seurat object for downstream analysis
seu <- CreateSeuratObject(res_mat, project=opt$sampleName,min.cells = 3, min.features = 200)

#store percentage mitochondria as metadata
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-")

#QC plot of features, counts and %mitochondria
v <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) +
  cowplot::theme_cowplot(font_size = 20)

g <- cowplot::plot_grid(barcodePlot,cellPlot,v,ncol = 1)

cowplot::save_plot(g, filename = opt$QCPlot, base_width=7, base_height=12,bg="white")

#filter the Seurat object to retain good quality cells
#TODO add more flexibility to options
seu <- subset(x = seu, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt >  -Inf & percent.mt < 5 )

#add the metaData
condition <- metaDataTable[metaDataTable$Sample==opt$sampleName,"Condition"]
seu$Condition <- condition

#save the seurat object for later use
saveRDS(seu,opt$seuratOutput)


