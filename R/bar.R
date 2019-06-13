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

bar <- function(Data.total.bar){
  par(mfrow = c(1, 1),mar = c(12,4,2,2))
  barplot(Data.total.bar,ylab="Concentration",main = "TP53 Concentration in Samples",xaxt="n")
  labels <- colnames(Data.total.bar)
  # x axis with ticks but without labels
  axis(1, labels = FALSE)
  # Plot x labs at default x position
  text(x =  1.13*seq_along(labels), y = par("usr")[3] - 1, srt = 45, adj = 1,
       labels = labels, xpd = TRUE, cex = .6)
}
