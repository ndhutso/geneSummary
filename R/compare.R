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
  geneName <-as.character(geneName)

  #get expression levels of TP53
  idx <- match(geneName,Data.DBTRG$rownames.D2a.)
  Data.DBTRG.TP53 <- t(as.numeric(Data.DBTRG[idx,-1]))
  colnames(Data.DBTRG.TP53) <- names(Data.DBTRG)[-1]
  Data.U87.TP53 <- t(as.numeric(Data.U87[idx,-1]))
  colnames(Data.U87.TP53) <- names(Data.U87)[-1]

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

  #write a description of the kind of test used and if there's assumed variance or not
  kable(testTableF,align=rep('c', 5))%>%
    kable_styling(bootstrap_options = c("striped", "hover"),full_width = F, position = "center",font_size = 10)%>%
    column_spec(1, bold = T, border_right = T)%>%
    column_spec(4,border_right = T)%>%
    column_spec(7,border_right = T)%>%
    add_header_above(c(" " = 1, "DBTRG" = 3, "U87" = 3, "Total" = 4))
  #Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
}
