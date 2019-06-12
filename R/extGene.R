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

extGene <- function(){
  dat1 <- readLines("data/GSE43452_family.soft.gz")
  idx_start <- grep("table_begin", dat1)
  idx_end <- grep("table_end", dat1)
  dat1a <- c()
  for(i in 1:length(idx_start))
  {
    dat1a <- c(dat1a, dat1[(idx_start[i]+1):(idx_end[i]-1)])
  }
  D1L <- lapply(dat1a, function(x){
    unlist(strsplit(x, split="\t"))
  })
  D1a <- data.frame(do.call(rbind, D1L[2:length(D1L)]))
  colnames(D1a) <- D1L[[1]]
  return(D1a)
}
