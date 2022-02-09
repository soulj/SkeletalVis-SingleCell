withCallingHandlers({
  
  library(methods) # Because Rscript does not always do this
  
  options('useFancyQuotes' = FALSE)
  
  suppressPackageStartupMessages(library("optparse"))
  suppressPackageStartupMessages(library("RGalaxy"))
  
  
  option_list <- list()
  
  option_list$accessionNumber <- make_option('--accessionNumber', type='character')
  option_list$sampleTable <- make_option('--sampleTable', type='character')
  option_list$downloadedFiles <- make_option('--downloadedFiles', type='character')

  opt <- parse_args(OptionParser(option_list=option_list))
  
  
  
  
  getRNASeqExpressionData <- function (accessionNumber = GalaxyCharacterParam(), sampleTable = GalaxyInputFile(required = T,
                                                                                                               formatFilter = "tabular"), downloadedFiles = GalaxyOutput("files","tabular"))
  {
    library(curl)
    library(tidyverse)
    
    downloadFastqFile <- function(fastq,md5){
      
      fileName <- basename(fastq)
      print(sprintf("downloading %s",fileName))
      curl_fetch_disk(fastq,fileName)
      #check the md5sum is correct - if not try downloading again
      if(tools::md5sum(fileName)!=md5) {
        print(sprintf("retrying %s",fileName))
        curl_fetch_disk(fastq,fileName)
        #if still wrong inform the user
        if(tools::md5sum(fileName)!=md5) {
          stop(sprintf("Error: MD5 sum of %s incorrect",fileName))
        }
      }
    }
    
    sampleTable <- read.delim(sampleTable)

    url <- paste0("http://www.ebi.ac.uk/ena/portal/api/filereport?accession=",accessionNumber,"&result=read_run&fields=study_accession,run_accession,fastq_ftp,fastq_md5&format=tsv")

    results <- read.delim(url,stringsAsFactors = FALSE)
    
    results <- results[results$fastq_ftp!= "", ] %>% separate_rows(fastq_ftp,fastq_md5,sep=";") %>% as.data.frame()
    results <- results[grep(paste(sampleTable$File, collapse = "|"),
                            results$run_accession),]
    results$fastq_ftp <- paste0("ftp://", results$fastq_ftp)
      
    mapply(downloadFastqFile,results$fastq_ftp,results$fastq_md5)

    }
    

  params <- list()
  for(param in names(opt))
  {
    if (!param == "help")
      params[param] <- opt[param]
  }
  
  setClass("GalaxyRemoteError", contains="character")
  wrappedFunction <- function(f)
  {
    tryCatch(do.call(f, params),
             error=function(e) new("GalaxyRemoteError", conditionMessage(e)))
  }
  
  
  suppressPackageStartupMessages(library(RGalaxy))
  do.call(getRNASeqExpressionData, params)
  
  ## end warning handler
}, warning = function(w) {
  cat(paste("Warning:", conditionMessage(w), "\n"))
  invokeRestart("muffleWarning")
})

    
