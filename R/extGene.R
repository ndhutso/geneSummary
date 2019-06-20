#'Extract Gene Annotation Matrix
#'
#'@description Extracts the gene annotation matrix
#'
#'@usage extSample()
#'
#'@author Nicholas Hutson
#'
#'@examples extSample()
#'
#'@export

extGene = function(data, geneSymbol = NA){

  name <- names(data)
  name <- str_remove(name, "_series_matrix.txt.gz")
  if(length(data)>1){
    sampleNote <- list()
    data3 <- list()
    for(n in 1:length(data)){
        data3[[n]] <- data[[n]]
        sampleNote[[n]] <- data.frame()
        sampleNote[[n]] <- rbind(sampleNote[[n]],pData(featureData(data3[[n]])))
        idxSym <- grep("Symbol", colnames(sampleNote[[n]]))
        if(!is.na(geneSymbol)){
          sampleNote[[n]] <- sampleNote[[n]][which(sampleNote[[n]][,idxSym]==geneSymbol,arr.ind = TRUE),]
        }
    }
  }else{
    data3 <- data
    sampleNote <- pData(featureData(data3[[1]]))
    idxSym <- grep("Symbol", colnames(sampleNote))
    if(!is.na(geneSymbol)){
      sampleNote <- sampleNote[which(sampleNote[,idxSym]==geneSymbol,arr.ind = TRUE),]
    }
  }
  return(list(name,sampleNote))
}
