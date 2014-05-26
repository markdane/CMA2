#Perform a rank product analysis on the scored topTable files
#in the lower directories
# Author Mark Dane

#Version .1
library("RankProd")

reduceToInts<-function(x,col="GeneID"){
  #Eliminate the unneeded columns
  
  x<-x[,c(col,"score","GeneID")]
  #Average the scores of any elements with more than one entry
  for(element in unique(x[,col])){
    x$score[x[,col]==element]<-mean(x$score[x[,col]==element])
  }
  #Reduce the dataframe to only one entry per col element
  x<-unique(x)
}

rankProd<-function(x){
  rankP<-1
  for(rank in x){
    rankP<-rankP*rank
  }
  return(rankP)
}

calcRP<-function(gene,ranks, col="GeneID"){
  temp<-sapply(ranks,rankings, gene=gene, col=col)
  rankProd<-rankProd(temp)
  return(rankProd)
}

rankings<-function(df,gene, col="GeneID"){
  which(df[,col] %in% gene)
}

#TODO: Get the p vaules for the rank products, not sure this is working yet
calcPs<-function(r,genes){
  #P value is based on the number of k experiments and n number of genes
  count<-length(genes %in% r)
  righttailgamma(r=rankProds,k=count*length(cellTypes),n=length(rankProds))
}

#------Begin Main----------
source("/Users/markdane/Documents/Thesis Process/RankProd/RankProdExact.R")
col="siRNAID"

#Get file names in working directory and sub directories that have the pattern in them
files<-dir(pattern = "median", full.names = TRUE, ignore.case = TRUE)

#Read the topTables into the a list of dataframes
ranks<-lapply(paste0(files,"/in/topTable.txt"),read.delim, stringsAsFactors=FALSE)

#Reduce the dataframes to the gene name and score columns,
#average technical replicates
ranks<-lapply(ranks,reduce,col=col)

#Calculate the rank products of the genes across all of the dataframes
rankProds<-sapply(ranks[[1]][[col]],calcRP,ranks=ranks, col=col)

#Get the p values for each rank product
pValues<-righttailgamma(r=rankProds,k=length(ranks),n=length(rankProds))

#Sort the genes by their rank products
pValues<-sort(pValues)

#Replace the siRNA names with the names of their target genes
cat(ranks[[1]]$GeneID[ranks[[1]]$siRNAID %in% names(pValues[pValues<=0.05])])
#TODO: test and troubleshoot
