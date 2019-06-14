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

bar <- function(data, geneSymbol=NA, dName=NA){

  #start by writing code assuming all variables are defined
  D2a <- extExp(data,geneSymbol,dName) #gene expression
  D2b <- extSample(data,dName) #sample annotations
  D1a <- extGene(data,dName) #gene annotations

  geneName <- D1a$ID[[match(geneSymbol,D1a$Symbol)]] #test with geneSymbol being a vector or add if statement to method that works with a vector
  Data.symbol <- data.frame(D2a[match(geneName,rownames(D2a)),])
  colnames(Data.symbol) <- as.character(D2b$title)

  Data.bar <- Data.symbol

  Data.bar <- tibble(Sample = colnames(Data.bar), Concentration = as.numeric(Data.bar))
  title <- paste(geneSymbol,"Concentration in Sampes") #going to have to change when geneSymbol is not defined?
  Data.bar %>%
    ggplot(aes(x = Sample, y=Concentration)) +
    geom_bar(stat="identity", position ="dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)) +
    ggtitle(title)
}





