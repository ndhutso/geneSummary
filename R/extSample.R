#'Extract Sample Matrix
#'
#'@description Extracts the sample annotation matrix
#'
#'@usage extSample(data,sampleName)
#'
#'@author Nicholas Hutson
#'
#'@examples data <- getGEO("GSE443452")
#'sampleName <- "DBTRG untreated rep1"
#'extSample(data,sampleName)
#'
#'@param data data imported using GEOquery
#'@param sampleName the sample name
#'
#'@export

extSample = function(data, sampleName = NA){

  name <- names(data)
  name <- str_remove(name, "_series_matrix.txt.gz")
  if(length(data)>1){
    sampleNote <- list()
    data3 <- list()
    for(n in 1:length(data)){
        data3[[n]] <- data[[n]]
        sampleNote[[n]] <- data.frame()
        sampleNote[[n]] <- rbind(sampleNote[[n]],pData(phenoData(data3[[n]]))) #might already take into account multiple data sets
        if(!is.na(sampleName)){
          sampleNote[[n]] <- sampleNote[[n]][match(sampleName,sampleNote[[n]]$title),]
        }
    }
  }else{
      data3 <- data
      sampleNote <- pData(phenoData(data3[[1]]))
      if(!is.na(sampleName)){
        sampleNote <- sampleNote[match(sampleName,sampleNote$title),]
      }
  }
  return(list(name,sampleNote))
}
