#'Run Shiny app
#'
#'@description Runs a Shiny app that allows the user to search the NCBI GEO database for data sets and creates data tables from them
#'
#'@usage runShiny()
#'
#'@examples runShiny()
#'
#'@author Nicholas Hutson
#'
#'@export
runShiny <- function(){

  appDir <- system.file("R", "Shiny.R", package = "geneSummary")

  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `geneSummary`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")

}
