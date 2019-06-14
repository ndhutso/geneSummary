data <- import()

D2a <- extExp(data) #gene expression

D2b <- extSample(data) #sample annotations

D1a <- extGene(data) #gene annotations

library(gsubfn)
list[Data.total.bar,Data.total.box,testTableF,testTibbleF3] <- analyze(D2a,D2b,D1a) #list[] not actually used in geneSummary package, just this script

bar(D1a,D2a,D2b)
box(D1a,D2a)

compare(D1a,D2a)
hist(D1a,D2a)

