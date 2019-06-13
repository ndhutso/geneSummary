data <- import()

D2a <- extExp(data)

D2b <- extSample(data)

D1a <- extGene(data)

library(gsubfn)
list[Data.total.bar,Data.total.box,testTableF,testTibbleF3] <- analyze(D2a,D2b,D1a) #list[] not actually used in geneSummary package, just this script

bar(Data.total.bar)
box(Data.total.box)

compare(testTableF)
hist(testTibbleF3)

