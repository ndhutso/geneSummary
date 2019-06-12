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
  list[D2c,dat2a] <- extSample()
  D1a <- extGene()
  list[Data.total.bar,Data.total.box,testTableF,testTibbleF3] <- analyze(D2a,D2c,dat2a,D1a)
  visualize(Data.total.bar,Data.total.box)
  compare(testTableF,testTibbleF3)
}
