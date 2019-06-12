#'Compare Data
#'
#'@description Compares the data
#'
#'@usage compare()
#'
#'@author Nicholas Hutson
#'
#'@examples compare()
#'
#'@export

compare <- function(){
  #write a description of the kind of test used and if there's assumed variance or not
  kable(testTableF,align=rep('c', 5))%>%
    kable_styling(bootstrap_options = c("striped", "hover"),full_width = F, position = "center",font_size = 10)%>%
    column_spec(1, bold = T, border_right = T)%>%
    column_spec(4,border_right = T)%>%
    column_spec(7,border_right = T)%>%
    add_header_above(c(" " = 1, "DBTRG" = 3, "U87" = 3, "Total" = 4))
  #Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

  #graph
  testTibbleF3 %>%
    ggplot(aes(x = Symbol, y=Mean, fill=Group.Mean)) +
    geom_bar(stat="identity", position ="dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                  position=position_dodge(.9)) +
    ggtitle("Average Gene Expression in Groups")
}
