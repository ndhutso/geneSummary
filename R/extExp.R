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
#'@import "Biobase"
#'@import "GEOquery"
#'@import "ggplot2"
#'@import "tibble"
#'@import "shiny"
#'
#'@export

extExp = function(data, geneSymbol=NA) {

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

      if(is.na(geneSymbol))
      {
          expData[[i]] <- as.data.frame(data2[[i]])
          geneSymbol <- data1[[i]][,idxSym]
          expData[[i]] <- add_column(expData[[i]],Symbol = geneSymbol, .before = colnames(expData[[i]])[[1]]) #somehow changing into a list here

      }else{

        geneName <- data1[[i]]$ID[match(geneSymbol,data1[[i]][,idxSym])] #might not account for multiple genes with same symbol

        if(length(geneName)==1){
          expData[[i]] <- data.frame(t(data2[[i]][match(geneName,rownames(data2[[i]])),]))
          expData[[i]] <- add_column(expData[[i]],Symbol = replicate(length(rownames(expData[[i]])), geneSymbol), .before = colnames(expData[[i]])[[1]])
        }else{
          geneSymbol <- data1[[i]][match(geneName,data1[[i]]$ID),idxSym] #have to set geneSymbol to length of columns bc some gene symbols repeat
          expData[[i]] <- data.frame(data2[[i]][match(geneName,rownames(data2[[i]])),]) #may have to use if statement for t() if there is one or more appearances of a symbol
          expData[[i]] <- add_column(expData[[i]],Symbol = geneSymbol, .before = colnames(expData[[i]])[[1]])
        }
      }
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

    if(is.na(geneSymbol))
    {
        expData <- as.data.frame(data2)
        geneSymbol <- data1[,idxSym]
        expData <- add_column(expData,Symbol = geneSymbol, .before = colnames(expData)[[1]])
    }else{

      geneName <- data1$ID[match(geneSymbol,data1[,idxSym])] #might not account for multiple genes with same symbol

      if(length(geneName)==1){
        expData <- data.frame(data2[match(geneName,rownames(data2)),])
        expData <- add_column(expData,Symbol = replicate(length(rownames(expData)), geneSymbol), .before = colnames(expData)[[1]])
      }else{
        geneSymbol <- data1[match(geneName,data1$ID),idxSym] #have to set geneSymbol to length of columns bc some gene symbols repeat
        expData <- data.frame(data2[match(geneName,rownames(data2)),]) #may have to use if statement for t() if there is one or more appearances of a symbol
        expData <- add_column(expData,Symbol = geneSymbol, .before = colnames(expData)[[1]])
      }
    }
  }
  return(list(name,expData))
}

