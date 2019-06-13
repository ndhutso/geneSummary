data <- import()

D2a <- extExp(data)

D2b <- extSample(data)

D1a <- extGene(data)

list[Data.total.bar,Data.total.box,testTableF,testTibbleF3] <- analyze(D2a,D2b,D1a) #have to have gsubfn loaded and environment clear, the package isn't loading gsubfn properly

res <- analyze(D2a,D2b,D1a)
names(res) <- c("Data.total.bar", "Data.total.box", "testTableF", "testTibbleF3")

with(res, dim(box))

with(res, visualize(Data.total.bar,Data.total.box))

with(res, compare(testTableF,testTibbleF3))

