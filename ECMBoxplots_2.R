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
  reps<-cbind(reps,tt[,c("ECM","GrowthFactor")])
  #Define a function that will plot all channels of a subset type
  boxFun<-function(x,reps=reps,type){
    simpleName<-gsub(pattern="Input Net, pin none, plate median",replacement="",x=experiment)
    #out<-boxplot(reps[,x]~reps[,type],main=paste(experiment,"\nChannel",x,type),cex.axis=.5,notch=TRUE)
    out<-boxplot(reps[,x]~reps[,type],main=paste(simpleName,"channel",x,"\n",type),cex.axis=.7,notch=TRUE,las=2,
                 yaxp=c(-10,10,4),ylim=c(-10,10))
    abline(h=0, col="blue")}
  #Call the function for different types of subsets
  for(type in c("ECM","GrowthFactor")){
    out<-sapply(1:nrChannels,boxFun,reps=reps,type=type)
  }
}#end ttBoxesPerSpot

ttBoxesPerCell<-function(dir,experiment,nrChannels=3){
  #read in the topTable from an ECM and GF analysis and make
  #boxplots of of all replicates by ECM and GF normed by the cellmask channel
  
  #TODO: Delete this line after debug
#      dir<- "/Users/markdane/Documents/Sophia MEMa 2/"
#      experiment<- "0 GY, Cells 10A, Input Net, pin none, plate median"
#      nrChannels=3
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
    out<-boxplot(reps[,x]~reps[,type],main=paste(simpleName,x,"\n",type),cex.axis=.7,notch=TRUE,las=2,
                 yaxp=c(-10,10,4),ylim=c(-10,10))
    abline(h=0, col="blue")
    medians<-out$stats[3,]
    #print(medians)
    #names(medians)<-names(reps[,type])
    #str(out$stats[3,])
    write.table(medians,file="temp.txt")
  return(out)
  }#End boxFunNorm
  
#----------Main call of ttBoxesPerCell---------------
  for(type in c("ECM","GrowthFactor")){
    out<-sapply(c("ch1Normch3", "ch2Normch3"),boxFunNorm,reps=reps,type=type)
  }
  
  return(out)
}#end ttBoxesPerCell


#----------Main-------------

setwd("/Users/markdane/Documents/Sophia MEMa 2/")
pdf(file="10AbyECMandGF.PDF")

spotBoxes<-sapply(c("0 GY, Cells 10A, Input Net, pin none, plate median",
                    "4 GY, Cells 10A, Input Net, pin none, plate median",
                    "8 GY, Cells 10A, Input Net, pin none, plate median"),FUN=ttBoxesPerSpot,dir=getwd())

#Plot boxplots grouped by ECM and Growth Factor of the per cell intensities to the pdf file.
cellBoxes<-sapply(c("0 GY, Cells 10A, Input Net, pin none, plate median",
         "4 GY, Cells 10A, Input Net, pin none, plate median",
         "8 GY, Cells 10A, Input Net, pin none, plate median"),FUN=ttBoxesPerCell,dir=getwd())
#Write a spreadsheet file of the median and standard deviations of the ECM and Growth Factor
#per cell intensities.

dev.off()

