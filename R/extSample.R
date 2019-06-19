#'Extract Sample Matrix
#'
#'@description Extracts the sample annotation matrix
#'
#'@usage extSample(data)
#'
#'@author Nicholas Hutson
#'
#'@examples extSample()
#'
#'@export

extSample = function(data, dName=NA){
  if(length(data)>1){
    sampleNote <- list()
    for(n in 1:length(data)){
      if(is.na(dName)){
        data3[[n]] <- data[[n]]
        sampleNote[[n]] <- data.frame()
        for(i in 1:length(data3[[n]]))
        {
          sampleNote[[n]] <- rbind(sampleNote[[n]],pData(phenoData(data3[[n]][[i]])))
        }
      }else if(length(dName)==1){
        data3[[n]] <- data[[grep(dName,names(data[[n]]),ignore.case = TRUE)]]
        sampleNote[[n]] <- pData(phenoData(data3[[n]]))
      }else{
        data3[[n]] <- data[[n]][[grep(dName,names(data[[n]]),ignore.case = TRUE)]]
        sampleNote[[n]] <- data.frame()
        for(i in 1:length(data3[[n]]))
        {
          sampleNote[[n]] <- rbind(sampleNote[[n]],pData(phenoData(data3[[n]][[i]])))
        }
      }
    }
  }else{
    if(is.na(dName)){
      data3 <- data
      sampleNote <- data.frame()
      for(i in 1:length(data3))
      {
        sampleNote <- rbind(sampleNote,pData(phenoData(data3[[i]])))
      }
    }else if(length(dName)==1){
      data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]]
      sampleNote <- pData(phenoData(data3))
    }else{
      data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]]
      sampleNote <- data.frame()
      for(i in 1:length(data3))
      {
        sampleNote <- rbind(sampleNote,pData(phenoData(data3[[i]])))
      }
    }
  }
  return(sampleNote)
}
