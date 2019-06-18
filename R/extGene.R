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
  return(sampleNote)
}
