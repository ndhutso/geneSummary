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
  #graphs histogram of first row of table expression data passed so one gene passed is preferred
  #browser()
  expData <- as.data.frame(t(tbl[,c(-1,-2)]))
  samples <- rownames(expData)
  expData <- add_column(expData,Sample = samples, .before = colnames(expData)[[1]])
  colnames(expData)[2] <- "Concentration"
  Data <- as_tibble(expData)

  graph <- Data %>%
    ggplot(aes(x = Sample, y=Concentration)) +
    geom_bar(stat="identity", position ="dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5))
  graph
  return(graph)
}
