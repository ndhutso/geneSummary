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
  analyze()
  visualize()
  compare()
}
