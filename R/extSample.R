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

extSample = function(data){

  name <- names(data)
  name <- str_remove(name, "_series_matrix.txt.gz")
  if(length(data)>1){
    sampleNote <- list()
    data3 <- list()
    for(n in 1:length(data)){
        data3[[n]] <- data[[n]]
        sampleNote[[n]] <- data.frame()
        sampleNote[[n]] <- rbind(sampleNote[[n]],pData(phenoData(data3[[n]]))) #might already take into account multiple data sets
    }
  }else{
      data3 <- data
      sampleNote <- pData(phenoData(data3[[1]]))
    }
  return(list(name,sampleNote))
}
