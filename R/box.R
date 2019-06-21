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


box <- function(D1a,D2a,D2b){

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

  ## Boxplot for the two groups, combine
  Data.total.box <- data.frame(c(replicate(4, "DBTRG"),replicate(8, "U87")),c(Data.total.bar[1,1:4],Data.total.bar[1,5:12]))
  colnames(Data.total.box) <- c("Group","TP53_Concentration")

  #graph
  par(mfrow = c(1, 1),mar=c(2,4,2,2))
  graph <- Data.total.box %>%
    ggplot(aes(Group, TP53_Concentration)) +
    geom_boxplot() +
    labs(title = "TP53 Concentration in Groups", y = "TP53 Concentration")
  graph
  return(graph)
}
