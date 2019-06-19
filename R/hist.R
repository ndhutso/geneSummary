#'Histogram of Data
#'
#'@description Compares the data in a histogram
#'
#'@usage compare()
#'
#'@author Nicholas Hutson
#'
#'@examples compare()
#'
#'@export

hist <- function(D1a,D2a){

  ## extract expression level for DBTRG group and U87 group separately - START TIDYING
  #get data of groups combined - change
  Data.DBTRG <- D2a[,1:5]
  Data.U87 <- data.frame(D2a[,1],D2a[,6:13])
  colnames(Data.U87)[1] = "Symbol"

  #list of gene names for TP53 exactly
  geneName <- D1a$ID[which(D1a$Symbol=="TP53",arr.ind=TRUE)]
  geneName <-as.character(geneName)

  #get expression levels of TP53
  idx <- match(geneName,rownames(Data.DBTRG))
  Data.DBTRG.TP53 <- t(as.numeric(Data.DBTRG[idx,-1]))
  colnames(Data.DBTRG.TP53) <- names(Data.DBTRG)[-1]
  Data.U87.TP53 <- t(as.numeric(Data.U87[idx,-1]))
  colnames(Data.U87.TP53) <- names(Data.U87)[-1]

  ## t test between the two groups
  test <- t.test(Data.DBTRG.TP53,Data.U87.TP53)

  #change sig to p-value
  test.Data.DBTRG.TP53 <- tibble(Mean.DBTRG = mean(Data.DBTRG.TP53),SD.DBTRG = sd(Data.DBTRG.TP53),Median.DBTRG = median(Data.DBTRG.TP53))
  test.Data.U87.TP53 <- tibble(Mean.U87 = mean(Data.U87.TP53),SD.U87 = sd(Data.U87.TP53),Median.U87 = median(Data.U87.TP53))

  total <- cbind(Data.DBTRG.TP53,Data.U87.TP53)
  test.total <- tibble(Mean = mean(total),SD = sd(total),Median = median(total))

  testTable <- cbind(ID_REF = geneName, test.Data.DBTRG.TP53,test.Data.U87.TP53,test.total,p.value = test$p.value[1])

  ##find top 10 variance and corresponding gene, add to summary table, visualize with ggplot
  geneName <- rownames(D2a)
  topVariance <- as_tibble(D2a) %>%
    mutate(Variance = rowVars(as.matrix(D2a[,-1]))) %>%
    mutate(ID_REF = geneName) %>%
    arrange(desc(Variance)) %>%
    slice(1:10)
  topVariance <- topVariance %>% select(ID_REF,colnames(topVariance))

  #separate data into groups, find mean median sd p-value, add
  #Data.U87 and Data.DBTRG are set
  geneName <- as.character(D1a$ID[match(topVariance$ID_REF,D1a$ID)])
  idx.D <- match(geneName,rownames(Data.DBTRG))
  idx.U <- match(geneName,rownames(Data.U87))
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
    mutate(Mean = rowMeans(topVariance[,3:14]),"SD" = rowSds(as.matrix(topVariance[,3:14])), Median = rowMedians(as.matrix(topVariance[,3:14]))) %>%
    mutate("p.value" = t(test[1,grep("p.value",colnames(test))]),ID_REF = topVariance$ID_REF)
  Data.total.Top <- Data.total.Top %>%
    select(-c(colnames(Data.total.Top[grep("GSM", colnames(Data.total.Top))]))) %>%
    select(-c(colnames(Data.total.Top[grep("Symbol", colnames(Data.total.Top))])))
  testTable2 <- plyr::join(Data.DBTRG.Top, Data.U87.Top, by = "ID_REF")
  testTable2 <- plyr::join(testTable2, Data.total.Top, by = "ID_REF")

  #create final table
  testTableF <- rbind(testTable,testTable2,stringsAsFactors = FALSE)

  #convert gene ID to actual symbols
  geneSymbol <- D1a$Symbol[match(testTableF$ID_REF, D1a$ID)] #assigns symbols of gene
  testTibbleF <- as_tibble(testTableF)
  colLabel <- c("Mean","SD", "Median")
  colLabel <- c(colLabel,colLabel,colLabel)
  colnames(testTableF)[1:10] <- c("Symbol",colLabel)

  #detach("package:plyr", unload = TRUE) #doesn't work with dplyr definition of rename
  #unloadNamespace("plyr")
  testTibbleF <- testTibbleF %>%
    dplyr::rename(Symbol=`ID_REF`) %>% #is breaking symbols
    mutate_if(is.character,as.numeric)
  testTibbleF$Symbol <- geneSymbol

  testTibbleF1 <- testTibbleF  %>% select(Symbol, Mean.DBTRG, Mean.U87)  %>%
    tidyr::gather(Group.Mean, Mean, -1) %>% #separates expression data based on group
    mutate(group=c(replicate(11,"DBTRG"),replicate(11,"U87")))

  testTibbleF2 <- testTibbleF  %>% select(Symbol, SD.DBTRG, SD.U87) %>%
    tidyr::gather(Group.SD, SD, -1) %>%
    mutate(group=c(replicate(11,"DBTRG"),replicate(11,"U87")))

  testTibbleF3 <- plyr::join(testTibbleF1,testTibbleF2, by=c("Symbol", "group"))

  #graph
  graph <- testTibbleF3 %>%
    ggplot(aes(x = Symbol, y=Mean, fill=Group.Mean)) +
    geom_bar(stat="identity", position ="dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2, position=position_dodge(.9)) +
    ggtitle("Average Gene Expression in Groups")
  graph
  return(graph)
}
