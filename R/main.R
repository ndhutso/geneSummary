#'Main
#'
#'@usage main()
#'
#'@author Nicholas Hutson
#'
#'@examples main()
#'
#'@export

main <- function(){
  import()
  D2a <- extExp()
  D2c <- extSample()
  D1a <- extGene()
  data <- analyze(D2a,D2c,D1a)
  visualize(data[[1]],data[[2]])
  compare(data[[3]],data[[4]])
}
