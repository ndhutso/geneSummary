library(shiny)
runExample("01_hello")

library(Biobase)
library(GEOquery)
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

##get sample annotations pData(phenoData(data[[1]]))
anotateSample = function(data, dName){
  if(is.na(dName)){
    data3 <- data
  }else if(length(dName)==1){
    data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]]
    sampleNote <- pData(phenoData(data3))
  }else{
    data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]]
    sampleNote <- data.frame()
    for(i in 1:length(data3))
    {
      sampleNote <- rbind(sampleNote,pData(phenoData(data3[[i]])))
    }
  }
  return(sampleNote)
}

##get gene annotations pData(featureData(data[[1]]))
annotateGene = function(data, dName){
  if(is.na(dName)){
    data3 <- data
  }else if(length(dName)==1){
    data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]]
    sampleNote <- pData(featureData(data3))
  }else{
    data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]]
    sampleNote <- data.frame()
    for(i in 1:length(data3))
    {
      sampleNote <- rbind(sampleNote,pData(featureData(data3[[i]])))
    }
  }
  return(geneNote)
}

##function to get gene expression data based on gene name - only searches first annotated data frame, could use grep to search names()
#might be able to merge geneExp and mergeExp
expGene = function(data, geneSymbol=NA, dName=NA) {

  if(is.na(dName)){
    data3 <- data
  }else if(length(dName)==1){
    data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]] #big list of all the data sets

    data1 <- pData(featureData(data3))
    data2 <- exprs(data3)
  }else{
    data3 <- data[[grep(dName,names(data),ignore.case = TRUE)]] #big list of all the data sets

    data1 <- data.frame()
    data2 <- data.frame()
    #have to use for loop or apply function to go through each data set, maybe use rbind to combine data
    for(i in 1:length(data3))
    {
      data1 <- rbind(data1, pData(featureData(data3[[i]])))
      data2 <- rbind(data2, exprs(data3[[i]]))
    }
  }

  if(is.na(geneSymbol))
  {
    geneName <- data1$ID
    geneSymbol <- data1$Symbol

    expData <- data.frame(data2[match(geneName,rownames(data2)),]) #may have to use if statement for t() if there is one or more appearances of a symbol
    expData <- add_column(expData,Symbol = replicate(length(rownames(expData)), geneSymbol), .before = colnames(expData)[[1]])
  } else {

    geneName <- data1$ID[match(geneSymbol,data1$Symbol)] #might not account for multiple genes with same symbol

    if(length(geneName)==1)
    {
      expData <- data.frame(t(data2[match(geneName,rownames(data2)),]))
      expData <- add_column(expData,Symbol = replicate(length(rownames(expData)), geneSymbol), .before = colnames(expData)[[1]])
    } else {
      geneSymbol <- data1$Symbol[match(geneName,data1$ID)] #have to set geneSymbol to length of columns bc some gene symbols repeat
      expData <- data.frame(data2[match(geneName,rownames(data2)),]) #may have to use if statement for t() if there is one or more appearances of a symbol
      expData <- add_column(expData,Symbol = geneSymbol, .before = colnames(expData)[[1]])
    }
  }

  return(expData)
}

data <- getGEO("GSE43452")
geneSymbol <- "TP53"
dName <- "GSE43452"
something(data,geneSymbol = c("TP53","STAR"),dName = "GSE43452")


##function to extract mean, std, median
geneSummary = function(expData)
{
  idx <- which(colnames(expData)!="Symbol")

  summary <- as_tibble(expData) %>%
    mutate_if(is.double,as.numeric) %>%
    mutate(Mean = rowMeans(expData[,idx])) %>%
    mutate(SD = rowSds(as.matrix(expData[,idx]))) %>%
    mutate(Median = rowMedians(as.matrix(expData[,idx])))

  return(summary)
}


geneExp = function(geneSymbol, data)
{
  data1 <- pData(featureData(data[[1]])) #D1a
  data2 <- exprs(data[[1]]) #D2a

  geneName <- data1$ID[which(data1$Symbol==geneSymbol,arr.ind = TRUE)]

  if(length(geneName)==1)
  {
    expData <- data.frame(t(data2[match(geneName,rownames(data2)),]))
    expData <- add_column(expData,Symbol = replicate(length(rownames(expData)), geneSymbol), .before = colnames(expData)[[1]])
  } else {
    expData <- data.frame(data2[match(geneName,rownames(data2)),]) #may have to use if statement for t() if there is one or more appearances of a symbol
    expData <- add_column(expData,Symbol = replicate(length(rownames(expData)), geneSymbol), .before = colnames(expData)[[1]])
    #turning single row data frames into lists
  }

  return(expData)
}


##test geneExp and geneSummary
phenoData(data[[1]]) #D2b
data <- getGEO("GSE43452")
geneSymbol <- "STAR"
expData <- geneExp(geneSymbol, data)
geneSummary(expData)

##function to merge data sets = combine

##function to merge and extract expression data
mergeExp = function(dataA, dataB, geneSymbol = NA)
{
  data3 <- combine(dataA, dataB)

  data1 <- data.frame()
  data2 <- data.frame()
  #have to use for loop or apply function to go through each data set, maybe use rbind to combine data
  for(i in 1:length(data3))
  {
    data1 <- rbind(data1, pData(featureData(data3[[i]])))
    data2 <- rbind(data2, exprs(data3[[i]]))
  }

  if(is.na(geneSymbol))
  {
    geneName <- data1$ID
    geneSymbol <- data1$Symbol

    expData <- data.frame(data2[match(geneName,rownames(data2)),]) #may have to use if statement for t() if there is one or more appearances of a symbol
    expData <- add_column(expData,Symbol = replicate(length(rownames(expData)), geneSymbol), .before = colnames(expData)[[1]])
  } else {

    geneName <- data1$ID[match(geneSymbol,data1$Symbol)] #might not account for multiple genes with same symbol

    if(length(geneName)==1)
    {
      expData <- data.frame(t(data2[match(geneName,rownames(data2)),]))
      expData <- add_column(expData,Symbol = replicate(length(rownames(expData)), geneSymbol), .before = colnames(expData)[[1]])
    } else {
      geneSymbol <- data1$Symbol[match(geneName,data1$ID)] #have to set geneSymbol to length of columns bc some gene symbols repeat
      expData <- data.frame(data2[match(geneName,rownames(data2)),]) #may have to use if statement for t() if there is one or more appearances of a symbol
      expData <- add_column(expData,Symbol = geneSymbol, .before = colnames(expData)[[1]])
    }
  }

  return(expData)
}

##test mergeExp
dataA <- getGEO("GSE43452")
dataB <- getGEO("GSE43452")
geneSymbol <- c("TP53","STAR", "DVL2","ECH1","TTN","CHFR")

mergeExp(dataA,dataB) #works, but vector takes up too much memory
mergeExp(dataA,dataB,geneSymbol)

##function to graph
#could do box plots comparing between the gene in the two different data sets




