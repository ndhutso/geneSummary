data <- getGEO("GSE43452")

D2a <- extExp(data) #gene expression

D2b <- extSample(data) #sample annotations

D1a <- extGene(data) #gene annotations

geneSymbol <- "TP53"
dName <- "GSE43452"

bar(data,"TP53","GSE43452")
box(data,"TP53","GSE43452")

compare(data,"TP53","GSE43452")
hist(data,"TP53","GSE43452")

