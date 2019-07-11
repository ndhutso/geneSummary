#'Filter Data Tables
#'
#'@description Filters data tables
#'
#'@usage filterTbl(tbl, type, long = FALSE, var, input)
#'
#'@author Nicholas Hutson
#'
#'@examples
#'
#'@export

filterTbl <- function(tbl, type, long = FALSE, var, input){ #OUTPUT ROW INDICES OR 0

  y <- match(var,colnames(tbl[[2]]))
  y <- y[!is.na(y)]
  w <- colnames(tbl[[2]])[y]

  #browser()

  x <- which(input!="",arr.ind = TRUE)
  y <- data.frame(tbl[[2]][,y])
  #get rid of all empty inputs
  y <- data.frame(y[,x])
  colnames(y) <- w[x]
  z <- input[x]

    if(type == "Gene Expression" & long == FALSE){
      #check if y input is not the first 2 column names of tbl
      #if so, search for rows where the sample gene expression fits the inputted logical expression
      #going to have to compare row indexes or logical vector with z vector
      if(identical(z,character(0))){
        z <- 0
      }else{

        #browser()

        a <- which(colnames(y) %in% colnames(tbl[[2]])[1:2], arr.ind = TRUE)#finds symbol and id inputs where there are inputs3

        if(!identical(a, integer(0))){
          y1 <- data.frame(y[,a])
          z1 <- z[a]

          z1 <- unlist(strsplit(z1,", ",fixed = TRUE)) #split up comma deliminated inputs
          z1 <- row.match(z1, y1) #could make this more generalized and complicated with grep

          #browser()

          y2 <- as.character(y[,-a])
          len <- dim(data.frame(y[,-a]))[2] #different way to find the number of elements
          z2 <- z[-a]
        }else{ #if there is no symbol or ID inputs (z1 = NULL or doesn't exist), what should z1 be
          y2 <- as.character(y)
          len <- dim(data.frame(y[,-a]))[2]
          z2 <- z
        }

        browser()

        if(!exists("z1")){
          if(len > 1){
            b <- list()
            for(i in 1:len){
              b[i] <- paste(y2[i], z2[i]) #should be like "y >4"
              b[i] <- parse(text = b[i])
            }
            z <- which(apply(data.frame(b),1,all)==TRUE, arr.ind = TRUE) #false when there's any false
          }else{
            b <- paste(y2, z2) #should be like "y >4"
            b <- parse(text = b)
            #need to have row indice output
            z <- which(unlist(lapply(b,eval))==TRUE, arr.ind = TRUE)
          }
        }else{
          if(identical(z2,character(0))){
            z <- z1
          }else if(len > 1){
            b <- list()
            for(i in 1:len){
              b[i] <- paste(y2[i], z2[i]) #should be like "y >4"
              b[i] <- parse(text = b[i])
            }
            z2 <- which(apply(data.frame(b),1,all)==TRUE, arr.ind = TRUE) #false when there's any false
            #now compare row indices
            if(z1 %in% z2){
              z <- z1
            }else{
              z <- 0
            }
          }else{
            b <- paste(y2, z2) #should be like "y >4"
            b <- parse(text = b)
            #need to have row indice output
            z2 <- which(unlist(lapply(b,eval))==TRUE, arr.ind = TRUE) #still have to parse with for or apply bc it's like putting y > 4 on an entire table
            #now compare row indices
            if(z1 %in% z2){
              z <- z1
            }else{
              z <- 0
            }
          }
        }
      } #should output the row indices for matching "Symbol" and "ID" inputs
    }else{
      if(identical(z,character(0))){
        z <- 0
      }else{
        z <- data.frame(strsplit(z,", ",fixed = TRUE)[[1]])
        z <- row.match(z, y) #could make this more generalized and complicated with grep
      }
    }
  return(z)
}
