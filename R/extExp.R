#'Extract Expression Matrix
#'
#'@description Extracts the gene expression matrix
#'
#'@usage extExp()
#'
#'@author Nicholas Hutson
#'
#'@examples extExp()
#'
#'@import "plyr"
#'@import "tidyverse"
#'@import "purrr"
#'@import "dplyr"
#'@import "sSeq"
#'@import "matrixStats"
#'@import "varhandle"
#'@import "knitr"
#'@import "kableExtra"
#'@import "magrittr"
#'@import "gsubfn"
#'@import "Biobase"
#'
#'@export

extExp <- function(){
  D2a <- read.csv("data/GSE43452_series_matrix.txt.gz", sep="\t", comment.char = "!")
  return(D2a)
}
