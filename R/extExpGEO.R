#'Extract Expression Matrix
#'
#'@description Extracts the gene expression matrix for GEO data
#'
#'@usage extExpGEO(data, geneSymbol=NA, long=FALSE)
#'
#'@author Nicholas Hutson
#'
#'@examples data <- getGEO("GSE443452")
#'geneSymbol <- "TP53"
#'long <- TRUE
#'extExp(data, geneSymbol, long)
#'
#'@param data data imported using GEOquery
#'@param geneSymbol the gene "name" or "symbol"
#'@param long a boolean value entered for the data returned to be in a long format
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
#'@import "ggplot2"
#'@import "tibble"
#'@import "shiny"
#'@import "stringr"
#'@import "DT"
#'@import "reshape2"
#'
#'@export

extExpGEO = function(data, geneSymbol=NA, long = FALSE) {

  name <- names(data)
  name <- str_remove(name, "_series_matrix.txt.gz")
  if(length(data)>1){
    data3 <- data

    data1 <- lapply(data3, fData)
    #CURSED CODE
    data2 <- lapply(data3, exprs)

    ##maybe put in for loop to go through every part of the list
    expData <- list()
    for(i in 1:length(data1)){

      idxSym <- grep("Symbol", colnames(data1[[i]]))
      idxName <- match("ID", colnames(data1[[i]]))

      if(is.na(geneSymbol))
      {
        expData[[i]] <- as.data.frame(data2[[i]])
        geneSymbol <- data1[[i]][,idxSym]
        geneName <- data1[[i]][,idxName]
        expData[[i]] <- add_column(expData[[i]],Symbol = geneSymbol, .before = colnames(expData[[i]])[[1]]) #somehow changing into a list here
        expData[[i]] <- add_column(expData[[i]],ID = geneName, .before = colnames(expData[[i]])[[1]])

      }else{

        geneName <- data1[[i]]$ID[match(geneSymbol,data1[[i]][,idxSym])] #might not account for multiple genes with same symbol

        if(length(geneName)==1){
          expData[[i]] <- data.frame(t(data2[[i]][match(geneName,rownames(data2[[i]])),]))
          expData[[i]] <- add_column(expData[[i]],Symbol = replicate(length(rownames(expData[[i]])), geneSymbol), .before = colnames(expData[[i]])[[1]])
          expData[[i]] <- add_column(expData[[i]],ID = replicate(length(rownames(expData[[i]])), geneName), .before = colnames(expData[[i]])[[1]])
        }else{
          geneSymbol <- data1[[i]][match(geneName,data1[[i]]$ID),idxSym] #have to set geneSymbol to length of columns bc some gene symbols repeat
          expData[[i]] <- data.frame(data2[[i]][match(geneName,rownames(data2[[i]])),]) #may have to use if statement for t() if there is one or more appearances of a symbol
          expData[[i]] <- add_column(expData[[i]],Symbol = geneSymbol, .before = colnames(expData[[i]])[[1]])
          expData[[i]] <- add_column(expData[[i]],ID = geneName, .before = colnames(expData[[i]])[[1]])
        }
      }
    }
    if(long){

    }
  }else{ #this section is good for one data set at a time
    data3 <- data

    data1 <- data.frame(lapply(data3, fData))
    #CURSED CODE
    data2 <- data.frame(lapply(data3, exprs))

    for(i in 1: length(colnames(data1))){
      colnames(data1)[i] <- strsplit(colnames(data1)[i],"gz.")[[1]][2]
    }

    for(i in 1: length(colnames(data2))){
      colnames(data2)[i] <- strsplit(colnames(data2)[i],"gz.")[[1]][2]
    }

    idxSym <- grep("Symbol", colnames(data1))
    idxName <- match("ID", colnames(data1))

    if(is.na(geneSymbol))
    {
      expData <- as.data.frame(data2)
      geneSymbol <- data1[,idxSym]
      geneName <- data1[,idxName]
      expData <- add_column(expData,Symbol = geneSymbol, .before = colnames(expData)[[1]])
      expData <- add_column(expData,ID = geneName, .before = colnames(expData)[[1]])
    }else{

      geneName <- data1$ID[match(geneSymbol,data1[,idxSym])] #might not account for multiple genes with same symbol

      if(length(geneName)==1){
        expData <- data.frame(data2[match(geneName,rownames(data2)),])
        expData <- add_column(expData,Symbol = replicate(length(rownames(expData)), geneSymbol), .before = colnames(expData)[[1]])
        expData <- add_column(expData,ID = replicate(length(rownames(expData)), geneName), .before = colnames(expData)[[1]])
      }else{
        geneSymbol <- data1[match(geneName,data1$ID),idxSym] #have to set geneSymbol to length of columns bc some gene symbols repeat
        expData <- data.frame(data2[match(geneName,rownames(data2)),]) #may have to use if statement for t() if there is one or more appearances of a symbol
        expData <- add_column(expData,Symbol = geneSymbol, .before = colnames(expData)[[1]])
        expData <- add_column(expData,ID = geneName, .before = colnames(expData)[[1]])
      }
    }
    if(long){
      #expression data is going to be a column so gene symbol and ID have to be repeated for the number of samples there are
      #first lets get symbol and ID vectors
      #samples are a repeat of column names for how many genes there are
      #exp data is going to be each row turn into a column and stacked
      expData <- melt(expData, variable.name = "Sample", value.name = "Expression")
    }
  }
  return(list(name,expData))
}
