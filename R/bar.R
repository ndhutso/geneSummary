#'Bar Plot Data
#'
#'@description Visualizes the data
#'
#'@usage bar()
#'
#'@author Nicholas Hutson
#'
#'@examples bar()
#'
#'@export

bar <- function(D1a,D2a,D2b){

  ## extract expression level for DBTRG group and U87 group separately - START TIDYING
  #get data of groups combined - change
  Data.DBTRG <- D2a[,1:5]
  Data.U87 <- data.frame(D2a[,1],D2a[,6:13])
  colnames(Data.U87)[1] = "Symbol"

  #list of gene names for TP53 exactly
  geneName <- D1a$ID[which(D1a$Symbol=="TP53",arr.ind=TRUE)]
  geneName <-as.character(geneName)

  #get expression levels of TP53
  idx <- match(geneName,rownames(Data.DBTRG))
  Data.DBTRG.TP53 <- t(as.numeric(Data.DBTRG[idx,-1]))
  colnames(Data.DBTRG.TP53) <- names(Data.DBTRG)[-1]
  Data.U87.TP53 <- t(as.numeric(Data.U87[idx,-1]))
  colnames(Data.U87.TP53) <- names(Data.U87)[-1]

  ## Barplot for TP53 in these samples, label samples
  Data.total.bar <- t(c(Data.DBTRG.TP53,Data.U87.TP53))
  colnames(Data.total.bar) <- as.character(D2b$title)

  #graph
  Data.total.bar <- tibble(Sample = colnames(Data.total.bar), Concentration = Data.total.bar[1,])
  graph <- Data.total.bar %>%
    ggplot(aes(x = Sample, y=Concentration)) +
    geom_bar(stat="identity", position ="dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)) +
    ggtitle("TP53 Concentration in Samples")
 graph
 return(graph)
}





