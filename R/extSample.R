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

extSample = function(data, dName="GSE43452"){
  if(is.na(dName)){
    data3 <- data
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
  return(sampleNote)
}
