#Perform a rank product analysis on the scored topTable files in the lower directories
# Author Mark Dane

#Version .1

reduce<-function(x){
  return(x[,c("GeneID","score")])
}

rankProd<-function(x){
  rankP<-1
  for(rank in x){
    rankP<-rankP*rank
  }
  return(rankP)
}

calcRP<-function(gene,ranks){
  temp<-sapply(ranks,rankings, gene=gene)
  rankProd<-rankProd(temp)
  return(rankProd)
}

rankings<-function(df,gene){
  which(df[,"GeneID"] %in% gene)
}

files<-dir(pattern = "(DAPI)", full.names = TRUE, ignore.case = TRUE)

ranks<-lapply(paste0(files,"/in/topTable.txt"),read.table, header=TRUE, stringsAsFactors=FALSE)

ranks<-lapply(ranks,reduce)

rankProds<-sapply(ranks[[1]][["GeneID"]],calcRP,ranks=ranks)

rankProds<-sort(rankProds)
