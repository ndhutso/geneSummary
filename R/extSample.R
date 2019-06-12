#'Extract Sample Matrix
#'
#'@description Extracts the sample annotation matrix
#'
#'@usage extSample()
#'
#'@author Nicholas Hutson
#'
#'@examples extSample()
#'
#'@export

extSample <- function(){
  dat2 <- readLines("data/GSE43452_series_matrix.txt.gz")
  idx_start <- grep("Sample", dat2)
  dat2a <- dat2[idx_start]
  D2c <- read.csv(textConnection(dat2a), sep="\t")
  D2c <- t(D2c)
  return(D2c)
}
