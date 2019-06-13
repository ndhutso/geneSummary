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


box <- function(Data.total.box)
{
  par(mfrow = c(1, 1),mar=c(2,4,2,2))
  boxplot(Data.total.box$DBTRG,Data.total.box$U87,names=c("DBTRG","U87"),ylab="Concentration",main = "TP53 Concentration in Groups")
}
