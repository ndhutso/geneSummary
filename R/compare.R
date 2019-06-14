#'Compare Data
#'
#'@description Compares the data in a table
#'
#'@usage compare()
#'
#'@author Nicholas Hutson
#'
#'@examples compare()
#'
#'@export

compare <- function(D1a,D2a){

  Data.DBTRG <- data.frame(rownames(D2a),D2a[,1:4])
  Data.U87 <- data.frame(rownames(D2a),D2a[,5:12])

  #list of gene names for TP53 exactly
  geneName <- D1a$ID[which(D1a$Symbol=="TP53",arr.ind=TRUE)]
  geneName <- as.character(geneName)

  #get expression levels of TP53
  idx <- match(geneName,Data.DBTRG$rownames.D2a.)
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
  idx.D <- match(geneName,Data.DBTRG$rownames.D2a.)
  idx.U <- match(geneName,Data.U87$rownames.D2a.)
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
  testTable2 <- plyr::join(Data.DBTRG.Top, Data.U87.Top, by = "ID_REF")
  testTable2 <- plyr::join(testTable2, Data.total.Top, by = "ID_REF")

  #create final table
  testTableF <- rbind(testTable,testTable2,stringsAsFactors = FALSE)

  #convert gene ID to actual symbols
  geneSymbol <- D1a$Symbol[match(testTableF$ID_REF, D1a$ID)] #assigns symbols of gene
  colLabel <- c("Mean","SD", "Median")
  colLabel <- c(colLabel,colLabel,colLabel)
  colnames(testTableF)[1:10] <- c("Symbol",colLabel)

  #write a description of the kind of test used and if there's assumed variance or not
  kable(testTableF,align=rep('c', 5))%>%
    kable_styling(bootstrap_options = c("striped", "hover"),full_width = F, position = "center",font_size = 10)%>%
    column_spec(1, bold = T, border_right = T)%>%
    column_spec(4,border_right = T)%>%
    column_spec(7,border_right = T)%>%
    add_header_above(c(" " = 1, "DBTRG" = 3, "U87" = 3, "Total" = 4))
  #Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
}
