options(Ncpus = 6)

#Install R libraries needed for pipeline
install.packages(c("curl", "httr"))
install.packages("devtools")

#install Bioconductor
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

#function to install libraries from bioconductor
installBioConductor <- function(libName) {
  if(libName %in% rownames(installed.packages()) == FALSE){
    BiocManager::install(libName,ask = FALSE)
  }}

#install packages defined in the text file
toinstall <- read.delim("install/libs.txt",header=F,stringsAsFactors = F)
sapply(toinstall[,1],installBioConductor)

#install liana from github
devtools::install_github('saezlab/liana',dependencies = FALSE)
