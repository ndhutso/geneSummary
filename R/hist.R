#'Histogram of Data
#'
#'@description Compares the data in a histogram
#'
#'@usage hist(D1a,D2a)
#'
#'@author Nicholas Hutson
#'
#'@examples data <- getGEO("GSE43452")
#'D2a <- extExp(data)[[2]]
#'D1a <- extGene(data)[[2]]
#'hist(D1a,D2a)
#'
#'@export

hist <- function(tbl){
  #graph
  browser()
  Data <- tibble(Sample = colnames(tbl)[c(-1,-2)], Concentration = tbl[,c(-1,-2)])
  graph <- Data %>%
    ggplot(aes(x = Sample, y=Concentration)) +
    geom_bar(stat="identity", position ="dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5))
  graph
  return(graph)
}
