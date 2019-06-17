data <- getGEO("GSE43452")

D2a <- extExp(data) #gene expression

D2b <- extSample(data) #sample annotations

D1a <- extGene(data) #gene annotations

bar(D1a,D2a,D2b)
box(D1a,D2a)

compare(D1a,D2a)
hist(D1a,D2a)

x <- "TP53, STAR"
geneSymbol <- strsplit(x,", ",fixed = TRUE)[[1]]
extExp(data,geneSymbol)
