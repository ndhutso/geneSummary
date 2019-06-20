data <- getGEO("GSE43452")

D2a <- extExp(data,"TP53")[[2]] #gene expression

D2b <- extSample(data,"DBTRG untreated rep1")[[2]] #sample annotations

D1a <- extGene(data,"TP53")[[2]] #gene annotations

bar(D1a,D2a,D2b)
box(D1a,D2a)

compare(D1a,D2a)
hist(D1a,D2a)

x <- "TP53"
geneSymbol <- strsplit(x,", ",fixed = TRUE)[[1]]
extExp(data, geneSymbol)

data <- getGEO("GSE40006") #list of 2, good for extExp test

D2a <- extExp(data)[[2]]

D2b <- extSample(data)[[2]] #sample annotations

D1a <- extGene(data)[[2]] #gene annotations
