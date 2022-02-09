library("optparse")

#Parse command line parameters
option_list <- list()

option_list$species <- make_option('--species', type='character', default = NULL,help="Species to make the reactome pathways for <human,mouse>")
option_list$reactomePathways <- make_option('--reactomePathways', type='character', default=NULL, help="File to save the pathways as")

opt <- parse_args(OptionParser(option_list=option_list))

library("reactome.db")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("tidyverse")


getReactomePathways <- function(species, outFile){
  
  #the reactome pathways are only provided with entrez IDs
  #so we need to convert them to gene symbols to match the data we have
  
  symbol2eg <- dplyr::case_when(
    species == "human" ~ "org.Hs.egSYMBOL",
    species == "mouse" ~ "org.Mm.egSYMBOL"
  )
  
  symbol2eg <- as.list(get(symbol2eg))
  
  #get eg to reactome pathway
  reactome2eg <- as.list(reactomePATHID2EXTID)
  
  speciesID <- dplyr::case_when(
    species == "human" ~ "R-HSA",
    species == "mouse" ~ "R-MMU"
  )
  
  #filter to obtain pathways for the selected species
  reactome2eg <- reactome2eg[grep(speciesID,names(reactome2eg))]
  
  #function to search through the pathway
  grepREACTOME <- function(id,mapkeys){
    unique(unlist(mapkeys[id],use.names=FALSE))
  }
  
  #convert the entrez ids to gene symbols for each pathway
  reactome <- lapply(reactome2eg,grepREACTOME,symbol2eg)
  
  #get the pathway names rather than the ids
  reactome2name <- as.list(reactomePATHID2NAME)
  reactomeNames <- sapply(names(reactome),grepREACTOME,reactome2name)
  names(reactome) <- reactomeNames
  
  #save the pathays for later use
  saveRDS(reactome,file=outFile)
  
}

getReactomePathways(opt$species,opt$reactomePathways)