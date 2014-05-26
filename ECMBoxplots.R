ttBoxesPerSpot<-function(dir,experiment,nrChannels=3){
  #read in the topTable from an ECM analysis and make
  #boxplots of of all replicates by ECM,GF and the combinations
  
  #TODO: Delete this line after debug
#    dir<- "/Users/markdane/Documents/Sophia MEMa 2/"
#    experiment<- "0 GY, Cells 10A, Input Net, pin none, plate median"
  par(mfrow=c(3,1))
  tt<-read.delim(paste0(dir,"/",experiment,"/in/topTable.txt"),stringsAsFactors=FALSE)
  names(tt)[names(tt)=="GeneID"]<-"ECM and GF"
  #Put all of the normalized intensities from the 4 replicates of each channel into a single column
  reps<-sapply(1:nrChannels,function(x){
    reps<-(c(tt[,paste0("normalized_r1_ch",x)],tt[,paste0("normalized_r2_ch",x)],
             tt[,paste0("normalized_r3_ch",x)],tt[,paste0("normalized_r4_ch",x)]))
  })
  
  #Add the ECM, Growth Factor and GeneID (Combined ECM and GF)annotations, repeating each set as needed
  reps<-cbind(reps,tt[,c("ECM","GrowthFactor","ECM and GF","Grid")])
  #Define a function that will plot all channels of a subset type
  boxFun<-function(x,reps=reps,type){
    out<-boxplot(reps[,x]~reps[,type],main=paste(experiment,"\nChannel",x,type),cex.axis=.5,notch=TRUE)
  abline(h=0, col="blue")}
  #Call the function for different types of subsets
  for(type in c("ECM","GrowthFactor")){
    out<-sapply(1:nrChannels,boxFun,reps=reps,type=type)
  }
}#end ttBoxesPerSpot

ttBoxesPerCell<-function(dir,experiment,nrChannels=3){
  #read in the topTable from an ECM analysis and make
  #boxplots of of all replicates by ECM,GF and the combinations
  
  #TODO: Delete this line after debug
  #    dir<- "/Users/markdane/Documents/Sophia MEMa 2/"
  #    experiment<- "0 GY, Cells 10A, Input Net, pin none, plate median"
  par(mfrow=c(2,1))
  tt<-read.delim(paste0(dir,"/",experiment,"/in/topTable.txt"),stringsAsFactors=FALSE)
  names(tt)[names(tt)=="GeneID"]<-"ECM and GF"
  #Put all of the normalized intensities from the 4 replicates of each channel into a single column
  reps<-sapply(1:nrChannels,function(x){
    reps<-(c(tt[,paste0("normalized_r1_ch",x)],tt[,paste0("normalized_r2_ch",x)],
             tt[,paste0("normalized_r3_ch",x)],tt[,paste0("normalized_r4_ch",x)]))
  })
  
  #Add the ECM, Growth Factor and GeneID (Combined ECM and GF)annotations, repeating each set as needed
  reps<-cbind(reps,tt[,c("ECM","GrowthFactor","ECM and GF","Grid")],
              ch1Normch3=reps[,1]-reps[,3],ch2Normch3=reps[,2]-reps[,3])
  
  boxFunNorm<-function(x,reps=reps,type){
    #Make box plots of the normed values grouped by type
    simpleName<-gsub(pattern="Input Net, pin none, plate median",replacement="",x=experiment)
    #     simpleName<-gsub(pattern="pin none,",replacement="",simpleName)
    out<-boxplot(reps[,x]~reps[,type],main=paste(simpleName,x,"\n",type),cex.axis=.29,notch=TRUE)
    abline(h=0, col="blue")
  return(out)
  }#End boxFunNorm
  
#----------Main---------------
  for(type in c("ECM","GrowthFactor")){
    out<-sapply(c("ch1Normch3", "ch2Normch3"),boxFunNorm,reps=reps,type=type)
  }
  
  return(out)
}#end ttBoxesPerCell

ttHitsPerCell<-function(dir,experiment,nrChannels=3){
  #read in the topTable from an ECM analysis and return
  #the hits of of all replicates grouped by ECM and GF
  
  medianByType<-function(type,subtype,data,cols){
    #Return the median of cols values that match type and subtype
    subtypes<-unique(data[,type])
    sapply(subtypes,function(subtype){
      select<-data[,type]==subtype
      sapply(data[select,cols],median)
    })
  }#medianByType
  
  madByType<-function(type,subtype,data,cols){
    #Return the median of cols values that match type and subtype
    subtypes<-unique(data[,type])
    sapply(subtypes,function(subtype){
      select<-data[,type]==subtype
      sapply(data[select,cols],mad)
    })
  }#madByType
  
  sdByType<-function(type,subtype,data,cols){
    #Return the median of cols values that match type and subtype
    subtypes<-unique(data[,type])
    sapply(subtypes,function(subtype){
      select<-data[,type]==subtype
      sapply(data[select,cols],sd)
    })
  }#sdByType
  
#   medianByType("GrowthFactor",reps,c("ch1Normch3","ch2Normch3"))
#   madByType("GrowthFactor",reps,c("ch1Normch3","ch2Normch3"))
#   sdByType("GrowthFactor",reps,c("ch1Normch3","ch2Normch3"))
  
  #TODO: Delete this line after debug
       dir<- "/Users/markdane/Documents/Sophia MEMa 2/"
       experiment<- "0 GY, Cells 10A, Input Net, pin none, plate median"
       nrChannels=3
  par(mfrow=c(2,1))
  tt<-read.delim(paste0(dir,"/",experiment,"/in/topTable.txt"),stringsAsFactors=FALSE)
  #Put all of the normalized intensities from the 4 replicates of each channel into a single column
  reps<-sapply(1:nrChannels,function(x){
    reps<-(c(tt[,paste0("normalized_r1_ch",x)],tt[,paste0("normalized_r2_ch",x)],
             tt[,paste0("normalized_r3_ch",x)],tt[,paste0("normalized_r4_ch",x)]))
  }#reps sapply function
               )
  
  #Add the ECM, Growth Factor annotations, repeating each set as needed
  #Add coluns that normalize the EdU by cell mask and CD44 by cellmask
  reps<-cbind(reps,tt[,c("ECM","GrowthFactor")],
              ch1Normch3=reps[,1]-reps[,3],ch2Normch3=reps[,2]-reps[,3])
  
  medians<-sapply(c("ECM","GrowthFactor"),medianByType,data=reps,
               subtype=subtype,cols=c("ch1Normch3","ch2Normch3"))
  
  mads<-sapply(c("ECM","GrowthFactor"),madByType,data=reps,
                  subtype=subtype,cols=c("ch1Normch3","ch2Normch3"))
  
  stdDevs<-sapply(c("ECM","GrowthFactor"),sdByType,data=reps,
                  subtype=subtype,cols=c("ch1Normch3","ch2Normch3"))

  return(temp=list(medians=medians,mads=mads,stdDevs=stdDevs))
}#ttHitsPerCell

     
#   hitsFunNorm<-function(x,reps=reps,type){
#     #Get the median and mad of the normed values grouped by type
#     simpleName<-gsub(pattern="Input Net, pin none, plate median",replacement="",x=experiment)
#     #     simpleName<-gsub(pattern="pin none,",replacement="",simpleName)
#     out<-boxplot(reps[,x]~reps[,type],main=paste(simpleName,x,"\n",type),cex.axis=.29,notch=TRUE)
#     sumDF<-summaryBy(formula=reps[,x]~reps[,type])
#   }
#   #   
#   #   #Call the function for per spot different types of subsets
#   #   for(type in c("ECM","GrowthFactor")){
#   #     out<-sapply(1:nrChannels,boxFunSpot,reps=reps,type=type)
#   #   }
#   #   
#   for(type in c("ECM","GrowthFactor")){
#     out<-sapply(c("ch1Normch3", "ch2Normch3"),hitsFunNorm,reps=reps,type=type)
#   }
#   return(sumDF)
#end ttHitsPerCell


#----------Main-------------

setwd("/Users/markdane/Documents/Sophia MEMa 2/")
pdf(file="10AbyECMandGFnormedbyCellMask.PDF")

# sapply(c("0 GY, Cells 10A, Input Net, pin none, plate median",
#                     "4 GY, Cells 10A, Input Net, pin none, plate median",
#                     "8 GY, Cells 10A, Input Net, pin none, plate median"),FUN=ttBoxesPerSpot,dir=getwd())

#Plot boxplots grouped by ECM and Growth Factor of the per cell intensities to the pdf file.
boxes<-sapply(c("0 GY, Cells 10A, Input Net, pin none, plate median",
         "4 GY, Cells 10A, Input Net, pin none, plate median",
         "8 GY, Cells 10A, Input Net, pin none, plate median"),FUN=ttBoxesPerCell,dir=getwd())
#Write a psreadsheet file of the median and standard deviations of the ECM and Growth Factor
#per cell intensities.
fileNames<-c("0 GY, Cells 10A, Input Net, pin none, plate median",
             "4 GY, Cells 10A, Input Net, pin none, plate median",
             "8 GY, Cells 10A, Input Net, pin none, plate median")
hits<-lapply(fileNames,FUN=ttHitsPerCell,dir=getwd())

#Create simpler names
simpleNames<-gsub(pattern="(Input Net, pin none, plate median)|[, ]",replacement="",x=fileNames)
simpleNames<-paste0("p",simpleNames)
names(hits)<-simpleNames
index<-1
for(radTreat in hits){
  ECMs<-sapply(radTreat,function(x){ x[[1]]
  })
  print(names(radTreats))
  print(ECMs)
}
  write.table(file=paste0(radTreats,"csv"),x=
  index<-index+1
  for(stat in radTreats){
    print(radTreats)
}
}

temp<-as.matrix(hits[[1]][[1]],hits[[1]][[2]])

write.csv(file=simpleNames[1],x=hits[[1]])
# hits<-sapply(c("0 GY, Cells 10A, Input Net, pin none, plate median",
#                "4 GY, Cells 10A, Input Net, pin none, plate median",
#                "8 GY, Cells 10A, Input Net, pin none, plate median"),FUN=ttHitsPerCell,dir=getwd())
hits<-sapply(c("0 GY, Cells 10A, Input Net, pin none, plate median"),FUN=ttHitsPerCell,dir=getwd(),USE.NAMES=TRUE)
dev.off()

