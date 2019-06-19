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

extGene = function(data, geneSymbol = NA, dName = NA){
  if(length(data)>1){
    sampleNote <- list()
    for(n in 1:length(data)){
      if(is.na(dName)){
        data3[[n]] <- data[[n]]
        sampleNote[[n]] <- data.frame()
        for(i in 1:length(data3[[n]]))
        {
          sampleNote[[n]] <- rbind(sampleNote[[n]],pData(featureData(data3[[n]][[i]])))
        }
        idxSym <- grep("Symbol", sampleNote[[n]])
        if(!is.na(geneSymbol)){
          sampleNote[[n]] <- sampleNote[[n]][which(sampleNote[[n]][,idxSym]==geneSymbol,arr.ind = TRUE),]
        }
      }else if(length(dName)==1){
        data3[[n]] <- data[[n]][[grep(dName,names(data[[n]]),ignore.case = TRUE)]]
        sampleNote[[n]] <- pData(featureData(data3[[n]]))
        idxSym <- grep("Symbol", sampleNote[[n]])
        if(!is.na(geneSymbol)){
          sampleNote[[n]] <- sampleNote[[n]][which(sampleNote[[n]][,idxSym]==geneSymbol,arr.ind = TRUE),]
        }
      }else{
        data3[[n]] <- data[[n]][[grep(dName,names(data[[n]]),ignore.case = TRUE)]]
        sampleNote[[n]] <- data.frame()
        for(i in 1:length(data3))
        {
          sampleNote[[n]] <- rbind(sampleNote[[n]],pData(featureData(data3[[n]][[i]])))
        }
        idxSym <- grep("Symbol", sampleNote[[n]])
        if(!is.na(geneSymbol)){
          sampleNote[[n]] <- sampleNote[[n]][which(sampleNote[[n]][,idxSym]==geneSymbol,arr.ind = TRUE),]
        }
      }
    }
  }else{
    if(is.na(dName)){
      data3 <- data
      sampleNote <- data.frame()
      for(i in 1:length(data3))
      {
        sampleNote <- rbind(sampleNote,pData(featureData(data3[[i]])))
      }
      idxSym <- grep("Symbol", sampleNote)
      if(!is.na(geneSymbol)){
        sampleNote <- sampleNote[which(sampleNote[,idxSym]==geneSymbol,arr.ind = TRUE),]
      }
    }else if(length(dName)==1){
      data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]]
      sampleNote <- pData(featureData(data3))
      idxSym <- grep("Symbol", sampleNote)
      if(!is.na(geneSymbol)){
        sampleNote <- sampleNote[which(sampleNote[,idxSym]==geneSymbol,arr.ind = TRUE),]
      }
    }else{
      data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]]
      sampleNote <- data.frame()
      for(i in 1:length(data3))
      {
        sampleNote <- rbind(sampleNote,pData(featureData(data3[[i]])))
      }
      idxSym <- grep("Symbol", sampleNote)
      if(!is.na(geneSymbol)){
        sampleNote <- sampleNote[which(sampleNote[,idxSym]==geneSymbol,arr.ind = TRUE),]
      }
    }
  }
  return(sampleNote)
}
