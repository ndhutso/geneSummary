% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extExpGEO.R
\name{extExpGEO}
\alias{extExpGEO}
\title{Extract Expression Matrix}
\usage{
extExpGEO(data, geneSymbol=NA, long=FALSE)
}
\arguments{
\item{data}{data imported using GEOquery}

\item{geneSymbol}{the gene "name" or "symbol"}

\item{long}{a boolean value entered for the data returned to be in a long format}
}
\description{
Extracts the gene expression matrix for GEO data
}
\examples{
data <- getGEO("GSE443452")
geneSymbol <- "TP53"
long <- TRUE
extExp(data, geneSymbol, long)

}
\author{
Nicholas Hutson
}
