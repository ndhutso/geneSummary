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

bar <- function(D1a,D2a){

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
  colnames(Data.total.bar) <- as.character(D2b$X.Sample_title)

  #graph
  par(mfrow = c(1, 1),mar = c(12,4,2,2))
  barplot(Data.total.bar,ylab="Concentration",main = "TP53 Concentration in Samples",xaxt="n")
  labels <- colnames(Data.total.bar)
  # x axis with ticks but without labels
  axis(1, labels = FALSE)
  # Plot x labs at default x position
  text(x =  1.13*seq_along(labels), y = par("usr")[3] - 1, srt = 45, adj = 1,
       labels = labels, xpd = TRUE, cex = .6)
}
