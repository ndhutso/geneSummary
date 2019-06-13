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

hist <- function(testTibbleF3){
  #graph
  testTibbleF3 %>%
    ggplot(aes(x = Symbol, y=Mean, fill=Group.Mean)) +
    geom_bar(stat="identity", position ="dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)) +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2, position=position_dodge(.9)) +
    ggtitle("Average Gene Expression in Groups")
}
