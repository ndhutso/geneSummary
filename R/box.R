#'Box Plot Data
#'
#'@description Visualizes the data
#'
#'@usage box(D1a,D2a)
#'
#'@author Nicholas Hutson
#'
#'@examples data <- getGEO("GSE43452")
#'D2a <- extExp(data)[[2]]
#'D1a <- extGene(data)[[2]]
#'box(D1a,D2a)
#'
#'@export


box <- function(tbl){
  #need to flesh out how to narrow down the options for this
  #browser()
  expData <- as.data.frame(t(tbl[,c(-1,-2)]))
  samples <- rownames(expData)
  expData <- add_column(expData,Sample = samples, .before = colnames(expData)[[1]])
  colnames(expData)[2] <- "Concentration"
  Data <- as_tibble(expData)

  par(mfrow = c(1, 1),mar=c(2,4,2,2))
  graph <- Data %>%
    ggplot(aes(Sample, Concentration)) +
    geom_boxplot()
  graph
  return(graph)
}
