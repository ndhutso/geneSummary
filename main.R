data <- import()

D2a <- extExp(data)

D2c <- extSample(data)

D1a <- extGene(data)

list[Data.total.bar,Data.total.box,testTableF,testTibbleF3] <- analyze(D2a,D2c,D1a)

visualize(Data.total.bar,Data.total.box)

compare(testTableF,testTibbleF3)

