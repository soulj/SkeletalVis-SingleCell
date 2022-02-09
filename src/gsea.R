#GSEA - find the pathways that are enriched in identified cell clusters
library("optparse")

#Parse command line parameters
option_list <- list()
option_list$diffExp <- make_option('--diffExp', type='character', default = NULL,help="TSV file with the differential expression")
option_list$reactomePathways <- make_option('--reactomePathways', type='character', default = NULL,help="Reactome pathways to use for enrichment")
option_list$enrichmentTable <- make_option('--enrichmentTable', type='character', default=NULL, help="Enrichment results")

opt <- parse_args(OptionParser(option_list=option_list))


library(fgsea)
library(ggplot2)
library(reactome.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(tidyverse)

#function to run the gene set enrichment
runGSEA <- function(diffExp,reactome) {
  #create a vector of log2FC and their gene names
  geneRanks <- diffExp$avg_log2FC
  names(geneRanks) <- diffExp$gene
  
  #do the gsea
  gseaRes <- fgseaMultilevel(reactome, geneRanks)

  #order the results
  gseaRes <- gseaRes[order(gseaRes$padj),]
  
  return(gseaRes)
  
}

#load the reactome pathways
reactome <- readRDS(opt$reactomePathways)

#read in the differential expression table
diffExp <- read.delim(opt$diffExp)

#run the GSEA on each cluster
gseaRes <- diffExp %>% group_by(cluster) %>%  group_map(~ runGSEA(.x,reactome)) %>%
  bind_rows(.id = "cluster") %>% filter(padj <=0.05) %>% as.data.frame()

#add the genes in the leading edge
gseaRes$leadingEdge <- sapply(gseaRes$leadingEdge,paste,collapse=" ")

#write out the results
write.table(gseaRes, file=opt$enrichmentTable, sep = "\t",row.names=F,quote = F,col.names = T)