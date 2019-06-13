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

compare <- function(testTableF){
  #write a description of the kind of test used and if there's assumed variance or not
  kable(testTableF,align=rep('c', 5))%>%
    kable_styling(bootstrap_options = c("striped", "hover"),full_width = F, position = "center",font_size = 10)%>%
    column_spec(1, bold = T, border_right = T)%>%
    column_spec(4,border_right = T)%>%
    column_spec(7,border_right = T)%>%
    add_header_above(c(" " = 1, "DBTRG" = 3, "U87" = 3, "Total" = 4))
  #Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
}
