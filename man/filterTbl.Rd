% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterTbl.R
\name{filterTbl}
\alias{filterTbl}
\title{Filter Data Tables}
\usage{
filterTbl(tbl, type, long = FALSE, var, input)
}
\description{
Returns the row indice of the corresponding values to the input.
}
\examples{
data <- getGEO("GSE443452")
geneSymbol <- "TP53"
long <- TRUE
tbl <- extExp(data, geneSymbol, long)
filterTbl(tbl, "Gene Expression", geneSymbol, long)

}
\author{
Nicholas Hutson
}
