#'Extract Expression Matrix
#'
#'@description Extracts the gene expression matrix
#'
#'@usage extExp(data, geneSymbol=NA, dName=NA)
#'
#'@author Nicholas Hutson
#'
#'@examples extExp()
#'
#'@import "plyr"
#'@import "tidyverse"
#'@import "purrr"
#'@import "dplyr"
#'@import "sSeq"
#'@import "matrixStats"
#'@import "varhandle"
#'@import "knitr"
#'@import "kableExtra"
#'@import "magrittr"
#'@import "gsubfn"
#'@import "Biobase"
#'@import "GEOquery"
#'
#'@export

extExp = function(data, geneSymbol=NA, dName=NA) {

  if(is.na(dName)){
    data3 <- data

    data1 <- lapply(data3, pData(featureData))
    data2 <- lapply(data3, exprs)
  }else if(length(dName)==1){
    data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]] #big list of all the data sets

    data1 <- pData(featureData(data3))
    data2 <- exprs(data3)
  }else{
    data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]] #big list of all the data sets

    data1 <- data.frame()
    data2 <- data.frame()
    #have to use for loop or apply function to go through each data set, maybe use rbind to combine data
    for(i in 1:length(data3))
    {
      data1 <- rbind(data1, pData(featureData(data3[[i]])))
      data2 <- rbind(data2, exprs(data3[[i]]))
    }
  }

  if(is.na(geneSymbol))
  {
    geneName <- data1$ID
    geneSymbol <- data1$Symbol

    expData <- data.frame(data2[match(geneName,rownames(data2)),]) #may have to use if statement for t() if there is one or more appearances of a symbol
    expData <- add_column(expData,Symbol = replicate(length(rownames(expData)), geneSymbol), .before = colnames(expData)[[1]])
  }else{

    geneName <- data1$ID[match(geneSymbol,data1$Symbol)] #might not account for multiple genes with same symbol

    if(length(geneName)==1)
    {
      expData <- data.frame(t(data2[match(geneName,rownames(data2)),]))
      expData <- add_column(expData,Symbol = replicate(length(rownames(expData)), geneSymbol), .before = colnames(expData)[[1]])
    }else{
      geneSymbol <- data1$Symbol[match(geneName,data1$ID)] #have to set geneSymbol to length of columns bc some gene symbols repeat
      expData <- data.frame(data2[match(geneName,rownames(data2)),]) #may have to use if statement for t() if there is one or more appearances of a symbol
      expData <- add_column(expData,Symbol = geneSymbol, .before = colnames(expData)[[1]])
    }
  }

  return(expData)
}
