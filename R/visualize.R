#'Visualize Data
#'
#'@description Visualizes the data
#'
#'@usage visualize()
#'
#'@author Nicholas Hutson
#'
#'@examples visualize()
#'
#'@export

visualize <- function(Data.total.bar,Data.total.box){
  par(mfrow = c(1, 1),mar = c(12,4,2,2))
  barplot(Data.total.bar,ylab="Concentration",main = "TP53 Concentration in Samples",xaxt="n")
  labels <- colnames(Data.total.bar)
  # x axis with ticks but without labels
  axis(1, labels = FALSE)
  # Plot x labs at default x position
  text(x =  1.13*seq_along(labels), y = par("usr")[3] - 1, srt = 45, adj = 1,
       labels = labels, xpd = TRUE, cex = .6)

  par(mfrow = c(1, 1),mar=c(2,4,2,2))
  boxplot(Data.total.box$DBTRG,Data.total.box$U87,names=c("DBTRG","U87"),ylab="Concentration",main = "TP53 Concentration in Groups")
}
