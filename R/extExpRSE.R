#'Extract Expression Matrix
#'
#'@description Extracts the gene expression matrix
#'
#'@usage extExp(data, geneSymbol=NA, long=FALSE)
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
#'@export

extExpRSE = function(data, geneSymbol=NA, long = FALSE) {

    data1 <- data.frame(rowData(data))
    #CURSED CODE

    data2 <- data.frame(assays(data)$counts)

    idxSym <- grep("symbol", colnames(data1))
    idxName <- grep("id", colnames(data1))

    if(is.na(geneSymbol))
    {
        expData <- data2
        geneSymbol <- as.character(data1[,idxSym])
        geneName <- as.character(data1[,idxName])
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
  return(expData)
}

