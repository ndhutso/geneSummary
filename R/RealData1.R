#take text files and organize them into matrix with column and row lables that have annotations
#samples as columns and genes as rows

#data1 <- read.csv('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE43nnn/GSE43452/soft/GSE43452_family.soft.gz')
#data2 <- read.csv('')
#data3 <- read.csv('')

#use "-" sign like this: data <- data[-c(1,2,3)] to remove rows

##OPTIMIZE USING BETTER METHODS, gene concentration table
#MOST OF THIS IS TRASH BUT USE AS EXAMPLE OF HOW SOME FUNCTIONS WORK
#data1 <- read.table('~/Downloads/GSE43452_family.soft.gz',fill=TRUE) #4 types of data for each gene, going to have to be 3D array 
#data1$V1 <- as.character(data1$V1)

#r <- NULL
#samples <- list()
#genes <- list()
#sBegin <- which(data1$V1=="!sample_table_begin", arr.ind = TRUE)
#sEnd <- which(data1$V1=="!sample_table_end", arr.ind = TRUE)

#for(i in 1:length(sBegin))
#{
#  sampledata <- data1[sBegin[i]:sEnd[i],]
#  sampledata <- sampledata[order(sampledata$V1),]
  
#  samples[[i]] <- sampledata[which(nchar(sampledata$V1)==12),2] 
#  names(samples)[i] <- paste("sample", i, sep = "")
  
#  r <- c(r, sBegin[i]:sEnd[i])
#}

#pBegin <- which(data1$V1=="!platform_table_begin", arr.ind = TRUE)
#pEnd <- which(data1$V1=="!platform_table_end", arr.ind = TRUE)

#platformTable <- data1[pBegin:pEnd,] #platform table with info on genes
#r <- c(r,pBegin:pEnd)

#data1a <- data1[-r,] #info on samples

#have to combine multiple data frames to one
#x <- eval(parse(text = paste("sample", p, sep = ""))) #how to associate chr to data.frame
#genes <- as.list(sampledata[which(nchar(sampledata$V1)==12),1])
#names(genes) <- sampledata[which(nchar(sampledata$V1)==12),1]

#dataFinal <- data.frame(lapply(samples, `[`))
#rownames(dataFinal) <- names(genes)

#now, add the extra information for the samples and genes
#each sample's data starts with ^SAMPLE and ends with !Sample_data_row_count
#sBegin <- which(data1a$V1=="^SAMPLE",arr.ind = TRUE)
#sEnd <- which(data1a$V1=="!Sample_data_row_count",arr.ind = TRUE)

#data2 <- read.table('~/Downloads/GSE43452_series_matrix.txt.gz',fill=TRUE)



##START OF QIANG DESTROYING MY WORK
#dat1 <- readLines("~/Downloads/GSE43452_family.soft.gz")
#idx_start <- grep("table_begin", dat1)
#idx_end <- grep("table_end", dat1)

#dat1a <- dat1[(idx_start[1]+1):(idx_end[1]-1)]
#head(dat1a,n=100)
# for loop = TRASH
#D1 <- c()
#for(i in 1:100){
#  r1 <- strsplit(dat1a[i], split = "\t")
#  D1 <- rbind(D1, unlist(r1))
#}
#colnames(D1) <- D1[1,]
#D1 <- D1[-1,]
# apply function
#D1L <- lapply(dat1a, function(x){ 
#  unlist(strsplit(x, split="\t"))  #splits large chr up by /t and then unlists to make a list
#})
#D1b <- do.call(rbind, D1L[2:length(D1L)])  #gets GENE TABLE, rbind binds rows of list to each other to create a table
#colnames(D1b) <- D1L[[1]] #first row of list is names so set column names to it

#D1L.30 <- D1L[lengths(D1L)==30]
#D1c <- do.call(rbind, D1L.30)

## 
#writeLines(dat1a, "/tmp/dat1a.txt")
#D1d <- read.table("/tmp/dat1a.txt", sep="\t")
#D1d <- read.csv("/tmp/dat1a.txt", sep="\t")

##
#D1e <- read.csv(textConnection(dat1a),sep="\t")



##START
library(tidyverse)
library(plyr)
library(purrr)
library(dplyr)
library(ggplot2)
library(sSeq) 
library(matrixStats)
library(varhandle)
library(knitr)
library(kableExtra)

##GENE TABLE
data1 <- read.table('~/Downloads/GSE43452_family.soft.gz',fill=TRUE) #can get data using GSEquery getGEO("GSE43452")
dat1 <- readLines("~/Downloads/GSE43452_family.soft.gz")

idx_start <- grep("table_begin", dat1)
idx_end <- grep("table_end", dat1)

dat1a <- c()
for(i in 1:length(idx_start))
{
  dat1a <- c(dat1a, dat1[(idx_start[i]+1):(idx_end[i]-1)])
}

D1L <- lapply(dat1a, function(x){ 
  unlist(strsplit(x, split="\t"))  #splits large chr up by /t and then unlists to make a list
})
D1a <- data.frame(do.call(rbind, D1L[2:length(D1L)]))  #gets GENE TABLE, rbind binds rows of list to each other to create a table
colnames(D1a) <- D1L[[1]] #first row of list is names so set column names to it



##OPTIMIZED TABLE FOR GENE CONCENTRATION
data2 <- read.table('~/Downloads/GSE43452_series_matrix.txt.gz',fill=TRUE) #onl to find what table certain data is on
dat2 <- readLines("~/Downloads/GSE43452_series_matrix.txt.gz") #turns data into giant chr

D2a <- read.csv("~/Downloads/GSE43452_series_matrix.txt.gz", sep="\t", comment.char = "!") #separates by "/t" and treats any line that starts with ! as a comment

dat2[grep("^!", dat2)] #uses a regular expression, ^, to find the lines that start with !



##SAMPLE TABLE - samples as rows and characteristics as columns
#find lines that start with "Sample"
idx_start <- grep("Sample", dat2)
idx_end <- grep("Sample", dat2)

dat2a <- c()

for(r in 1:length(idx_start))
{
  dat2a <- c(dat2a, dat2[(idx_start[r]+1):(idx_end[r]-1)])
}

dat2a <- dat2[idx_start]

D2L <- lapply(dat2a, function(x){ 
  unlist(strsplit(x, split="\t"))  #splits large chr up by /t and then unlists to make a list
})
D2b <- do.call(rbind, D2L[1:length(D2L)])  #gets GENE TABLE, rbind binds rows of list to each other to create a table
D2b <- t(D2b) #use reshape package to transpose rows and columns
colnames(D2b) <- D2b[1,] #first row of list is names so set column names to it
D2b <- data.frame(D2b[-1,])

D2c <- read.csv(textConnection(dat2a), sep="\t")
D2c <- t(D2c)



## extract expression level for DBTRG group and U87 group separately - START TIDYING
#get data of groups combined - change
Data.DBTRG <- data.frame(D2a[,1],D2a[,2:5]) 
Data.U87 <- data.frame(D2a[,1],D2a[,6:13])

#list of gene names for TP53 exactly
geneName <- D1a$ID[which(D1a$Symbol=="TP53",arr.ind=TRUE)]
geneName <-as.character(geneName) 

#get expression levels of TP53
idx.D <- match(geneName,Data.DBTRG$D2a...1.)
idx.U <- match(geneName,Data.U87$D2a...1.)

Data.DBTRG.TP53 <- t(as.numeric(Data.DBTRG[idx.D,-1])) 
colnames(Data.DBTRG.TP53) <- names(Data.DBTRG)[-1]

Data.U87.TP53 <- t(as.numeric(Data.U87[idx.U,-1]))
colnames(Data.U87.TP53) <- names(Data.U87)[-1]

## Barplot for TP53 in these samples, label samples
par(mfrow = c(1, 1),mar = c(12,2,2,2))
Data.total.bar <- t(c(Data.DBTRG.TP53,Data.U87.TP53))
colnames(Data.total.bar) <- as.character(D2b$X.Sample_title)
barplot(Data.total.bar,las=2,cex.names = .6)


## Boxplot for the two groups, combine
par(mfrow = c(1, 1))
Data.total.box <- data.frame(c(Data.total.bar[1,1:4],c(NA,NA,NA,NA)),Data.total.bar[1,5:12])
rownames(Data.total.box) <- NULL
colnames(Data.total.box) <- c("DBTRG","U87")
boxplot(Data.total.box$DBTRG,Data.total.box$U87,names=c("DBTRG","U87"))


## t test between the two groups
test <- t.test(Data.DBTRG.TP53,Data.U87.TP53)


## sumarized table: https://www.statalist.org/forums/forum/general-stata-discussion/general/1395253-descriptive-statistics-table-generation
#change sig to p-value
test.Data.DBTRG.TP53 <- c(mean(Data.DBTRG.TP53),sd(Data.DBTRG.TP53),median(Data.DBTRG.TP53)) 
test.Data.U87.TP53 <- c(mean(Data.U87.TP53),sd(Data.U87.TP53),median(Data.U87.TP53))
total <- c(Data.DBTRG.TP53,Data.U87.TP53)
test.total <- c(mean(total),sd(total),median(total)) 

testTable <- c(geneName, test.Data.DBTRG.TP53,test.Data.U87.TP53,test.total,test$p.value[1])
colLabel <- c("ID_REF","Mean.DBTRG","SD.DBTRG","Median.DBTRG","Mean.U87","SD.U87","Median.U87","Mean","SD","Median","p.value")

testTable <- t(testTable)
colnames(testTable) <- colLabel

#END TIDYING


##find top 10 variance and corresponding gene, add to summary table, visualize with ggplot
#should be able to use var function

D2a <- read.csv("~/Downloads/GSE43452_series_matrix.txt.gz", sep="\t", comment.char = "!")

topVariance <- as_tibble(D2a) %>%
  mutate(Variance = rowVars(as.matrix(D2a[,-1]))) %>%
  arrange(desc(Variance)) %>%
  slice(1:10)

#separate data into groups, find mean median sd p-value, add 
#Data.U87 and Data.DBTRG are set
geneName <- D1a$ID[match(topVariance$ID_REF,D1a$ID)]
geneName <-as.character(geneName) 

idx.D <- match(geneName,Data.DBTRG$D2a...1.)
idx.U <- match(geneName,Data.U87$D2a...1.)

Data.DBTRG.Top <- Data.DBTRG[idx.D,-1]
Data.U87.Top <- Data.U87[idx.U,-1]

test <- c()
for(i in 1:10) {
  test <- c(test, t.test(Data.DBTRG.Top[i,], Data.U87.Top[i,]))
}

test <- data.frame(test)

Data.DBTRG.Top <- as_tibble(Data.DBTRG.Top) %>%
  mutate(ID_REF = topVariance$ID_REF) %>%
  mutate("Mean.DBTRG" = rowMeans(Data.DBTRG.Top),"SD.DBTRG" = rowSds(as.matrix(Data.DBTRG.Top)), "Median.DBTRG" = rowMedians(as.matrix(Data.DBTRG.Top))) %>%
  select(-c(colnames(Data.DBTRG.Top[grep("GSM", colnames(Data.DBTRG.Top))])))

Data.U87.Top <- as_tibble(Data.U87.Top ) %>% 
  mutate(ID_REF = topVariance$ID_REF) %>%
  mutate("Mean.U87" = rowMeans(Data.U87.Top),"SD.U87" = rowSds(as.matrix(Data.U87.Top)), "Median.U87" = rowMedians(as.matrix(Data.U87.Top))) %>%
  select(-c(colnames(Data.U87.Top[grep("GSM", colnames(Data.U87.Top))])))


Data.total.Top <- as_tibble(topVariance) %>%
  select(2:13) %>%
  mutate(Mean = rowMeans(topVariance[,2:13]),"SD" = rowSds(as.matrix(topVariance[,2:13])), Median = rowMedians(as.matrix(topVariance[,2:13]))) %>%
  mutate("p.value" = t(test[1,grep("p.value",colnames(test))]),ID_REF = topVariance$ID_REF) 
Data.total.Top <- Data.total.Top %>%
  select(-c(colnames(Data.total.Top[grep("GSM", colnames(Data.total.Top))])))

testTable2 <- join(Data.DBTRG.Top, Data.U87.Top, by = "ID_REF") 
testTable2 <- join(testTable2, Data.total.Top, by = "ID_REF")

#create final table
testTableF <- rbind(testTable,testTable2,stringsAsFactors = FALSE)  

#convert gene ID to actual symbols
geneSymbol <- D1a$Symbol[match(testTableF$ID_REF, D1a$ID)]
testTableF$ID_REF <- geneSymbol
testTibbleF <- as_tibble(testTableF)
colLabel <- c("Mean","SD", "Median")
colLabel <- c(colLabel,colLabel,colLabel)
colnames(testTableF)[1:10] <- c("Symbol",colLabel)

#graph
detach("package:plyr") #doesn't work with dplyr definition of rename
testTibbleF <- testTibbleF %>% 
  rename(Symbol=`ID_REF`) %>%
  mutate_if(is.character,as.numeric)
library(plyr)
  
#testTibbleF$Symbol <- factor(testTibbleF$Symbol, testTibbleF$Symbol[order(testTibbleF$Mean)]) #changes factor level so they graph in the way the tibble is arranged

##add error bars to bar plot
testTibbleF1 <- testTibbleF  %>% select(Symbol, Mean.DBTRG, Mean.U87)  %>%
  gather(Group.Mean, Mean, -1) %>% #separates expression data based on group 
  mutate(group=c(replicate(11,"DBTRG"),replicate(11,"U87")))  

testTibbleF2 <- testTibbleF  %>% select(Symbol, SD.DBTRG, SD.U87) %>%
  gather(Group.SD, SD, -1) %>%
  mutate(group=c(replicate(11,"DBTRG"),replicate(11,"U87")))

testTibbleF3 <- join(testTibbleF1,testTibbleF2, by=c("Symbol", "group")) 

testTibbleF3 %>%
  ggplot(aes(x = Symbol, y=Mean, fill=Group.Mean)) + 
  geom_bar(stat="identity", position ="dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))