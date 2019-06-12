#'Import Data
#'
#'@description Imports the data
#'
#'@usage import()
#'
#'@author Nicholas Hutson
#'
#'@examples import()
#'
#'@export

import <- function(){
  dir.create("data/")
  download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE43nnn/GSE43452/soft/GSE43452_family.soft.gz", destfile = "data/GSE43452_family.soft.gz")
  download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE43nnn/GSE43452/matrix/GSE43452_series_matrix.txt.gz", destfile = "data/GSE43452_series_matrix.txt.gz")
}
