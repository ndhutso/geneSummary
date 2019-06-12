#'Analyze data
#'
#'@description Analyzes the data
#'
#'@usage analyze()
#'
#'@author Nicholas Hutson
#'
#'@examples analyze()
#'
#'@export

analyze <- function(){
  D2L <- lapply(dat2a, function(x){
    unlist(strsplit(x, split="\t"))  #splits large chr up by /t and then unlists to make a list
  })
  D2b <- do.call(rbind, D2L[1:length(D2L)])  #gets GENE TABLE, rbind binds rows of list to each other to create a table
  D2b <- t(D2b) #use reshape package to transpose rows and columns
  colnames(D2b) <- D2b[1,] #first row of list is names so set column names to it
  D2b <- data.frame(D2b[-1,])
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
  ## Boxplot for the two groups, combine
  par(mfrow = c(1, 1))
  Data.total.box <- data.frame(c(Data.total.bar[1,1:4],c(NA,NA,NA,NA)),Data.total.bar[1,5:12])
  rownames(Data.total.box) <- NULL
  colnames(Data.total.box) <- c("DBTRG","U87")
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
  detach("package:plyr") #doesn't work with dplyr definition of rename
  testTibbleF <- testTibbleF %>%
    rename(Symbol=`ID_REF`) %>%
    mutate_if(is.character,as.numeric)
  library(plyr)
  testTibbleF1 <- testTibbleF  %>% select(Symbol, Mean.DBTRG, Mean.U87)  %>%
    gather(Group.Mean, Mean, -1) %>% #separates expression data based on group
    mutate(group=c(replicate(11,"DBTRG"),replicate(11,"U87")))

  testTibbleF2 <- testTibbleF  %>% select(Symbol, SD.DBTRG, SD.U87) %>%
    gather(Group.SD, SD, -1) %>%
    mutate(group=c(replicate(11,"DBTRG"),replicate(11,"U87")))

  testTibbleF3 <- join(testTibbleF1,testTibbleF2, by=c("Symbol", "group"))
}
