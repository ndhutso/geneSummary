#'Box Plot Data
#'
#'@description Visualizes the data
#'
#'@usage box()
#'
#'@author Nicholas Hutson
#'
#'@examples box()
#'
#'@export


box <- function(D1a,D2a){

  ## extract expression level for DBTRG group and U87 group separately - START TIDYING
  #get data of groups combined - change
  Data.DBTRG <- data.frame(rownames(D2a),D2a[,1:4])
  Data.U87 <- data.frame(rownames(D2a),D2a[,5:12])

  #list of gene names for TP53 exactly
  geneName <- D1a$ID[which(D1a$Symbol=="TP53",arr.ind=TRUE)]
  geneName <-as.character(geneName)

  #get expression levels of TP53
  idx <- match(geneName,Data.DBTRG$rownames.D2a.)
  Data.DBTRG.TP53 <- t(as.numeric(Data.DBTRG[idx,-1]))
  colnames(Data.DBTRG.TP53) <- names(Data.DBTRG)[-1]
  Data.U87.TP53 <- t(as.numeric(Data.U87[idx,-1]))
  colnames(Data.U87.TP53) <- names(Data.U87)[-1]

  ## Barplot for TP53 in these samples, label samples
  Data.total.bar <- t(c(Data.DBTRG.TP53,Data.U87.TP53))
  colnames(Data.total.bar) <- as.character(D2b$title)

  ## Boxplot for the two groups, combine
  Data.total.box <- data.frame(c(Data.total.bar[1,1:4],c(NA,NA,NA,NA)),Data.total.bar[1,5:12])
  rownames(Data.total.box) <- NULL
  colnames(Data.total.box) <- c("DBTRG","U87")

  #graph
  par(mfrow = c(1, 1),mar=c(2,4,2,2))
  boxplot(Data.total.box$DBTRG,Data.total.box$U87,names=c("DBTRG","U87"),ylab="Concentration",main = "TP53 Concentration in Groups")
}
