#Functions that summarize technical replicates within a plate
#Author Mark Dane

repCVs<-function(dir,experiment,nrChannels=3){
  #read in the topTable from a MEMA analysis and add the mean, sd and cv across the 
  #biological and technical replicates
  #Do all calculations on the raw data assuming it has not been log transformed
  #TODO: add a paramter indicate logged raw values and backtransform if true
  
  #**********DEBUG**********: Comment out these line after debug
  #   dir<- "/Users/markdane/Documents/Sophia MEMa 2/"
  #   experiment<- "0 GY, Cells 10A, Input Net, pin none, plate median"
  #   par(mfrow=c(3,1))
  tt<-read.delim(paste0(dir,"/",experiment,"/in/topTable.txt"),stringsAsFactors=FALSE)
  names(tt)[names(tt)=="GeneID"]<-"ECMandGF"
  
  #There are 18 nonreplicated columns and 9 column types with replicated values
  nrReplicates<-(dim(tt)[2]-18)/9
  
  #Calculate the coefficient of variation for each biological replicate on a spot basis
  for(channel in 1:nrChannels){
    #create names for the columns
    bioRepMean<-paste0("bioRepMeanCh",channel)
    bioRepSD<-paste0("bioRepSDCh",channel)
    bioRepCV<-paste0("bioRepCVCh",channel)
    techRepMean<-paste0("techRepMeanCh",channel)
    techRepSD<-paste0("techRepSDCh",channel)
    techRepCV<-paste0("techRepCVCh",channel)
    
    #Calculate the mean of the normalized intensities of the 4 replicates of each spot and store them in a new column
    tt[,bioRepMean]<-apply(X=cbind(tt[,paste0("normalized_r1_ch",channel)],tt[,paste0("normalized_r2_ch",channel)],
                                   tt[,paste0("normalized_r3_ch",channel)],tt[,paste0("normalized_r4_ch",channel)]),
                           MARGIN=1,FUN=mean,na.rm=TRUE)
    
    #Calculate the sd of the normalized intensities of the 4 replicates of each spot and store them in a new column
    tt[,bioRepSD]<-apply(X=cbind(tt[,paste0("normalized_r1_ch",channel)],tt[,paste0("normalized_r2_ch",channel)],
                                 tt[,paste0("normalized_r3_ch",channel)],tt[,paste0("normalized_r4_ch",channel)]),
                         MARGIN=1,FUN=sd,na.rm=TRUE)
    
    #Calculate the CV of the normalized intensities of the 4 replicates of each spot and store them in a new column
    tt[,bioRepCV]<-tt[,bioRepSD]/abs(tt[,bioRepMean])
    
    #Calculate the technical replicate values within a plate.
    #Put the same value in the row of each replicate
    
    for(r in 1:nrReplicates){#Go through all the replicate plates in the screen
      for(name in unique(tt$ECMandGF)){#Get the name of one of the ECM+GF combinations
        techRepMean<-paste0("techRep",r,"MeanCh",channel)
        techRepSD<-paste0("techRep",r,"SDCh",channel)
        techRepCV<-paste0("techRep",r,"CVCh",channel)
        techRepCVCount<-paste0("techRep",r,"CVCountCh",channel)
        #Put the count of the non NA values of the current channel and current ECM+GF
        #into a new techRepXCVCountChX column
        tt[tt$ECMandGF==name,techRepCVCount]<-sum(!is.na(tt[tt$ECMandGF==name,paste0("normalized_r",r,"_ch",channel)]))
        
        #Put the mean of the raw values of the current channel,replicate and ECM+GF
        #into a new techRepXMeanchX column
        tt[tt$ECMandGF==name,techRepMean]<-mean(tt[tt$ECMandGF==name,
                                                   paste0("raw_r",r,"_ch",channel)],
                                                na.rm=TRUE)
        #Put the sd of the normalized values of the current channel and current ECM+GF
        #into a new techRepXSDchX column
        tt[tt$ECMandGF==name,techRepSD]<-sd(tt[tt$ECMandGF==name,
                                               paste0("raw_r",r,"_ch",channel)],
                                            na.rm=TRUE)
        #Put the cv of the normalized values of the current channel and current ECM+GF
        #into a new techRepXCVchX column
        tt[tt$ECMandGF==name,techRepCV]<-tt[tt$ECMandGF==name,techRepSD]/abs(tt[tt$ECMandGF==name,techRepMean])
      }#end for name in tt$ECMandGF
    }#End for r in nrReplicates
  }#End per channel calculations
  return(tt)
}#End repCVs function

# calcRatios()function(dir,experiment,nrChannels=3){
#   #read in the topTable from a MEMA analysis and add the ratios of the normalized values of 
#   #ch1/ch3 and ch2/ch3
#-----------Main-----------------------

experimentName<-"10A 8 Gy Analysis"
simpleName<-gsub(pattern="[, ]",x=experimentName,replacement="")

#Start a timer for optimization
ptb<-proc.time()

#Get a topTable and add the biological and technical replicate calculations
ttReps<-repCVs(dir= "/Users/markdane/Documents/Sophia MEMa 2/",
               experiment="8 GY, Cells 10A, Input Net, pin none, plate median", nrChannels=3)
#Write the topTable with replicate calculations to disk
write.table(x=ttReps,file="ttReps.txt",quote=FALSE,sep="\t",row.names=FALSE)

print(proc.time()-ptb)#takes 40-65 seconds on first implementation

ECMs<-unique(ttReps$ECM)
GFs<-unique(ttReps$GrowthFactor)
ECMandGFs<-unique(ttReps$ECMandGF)
ECMDF<-data.frame(ECM=ECMs)
GFDF<-data.frame(GF=GFs)
ECMandGFDF<-data.frame(ECMandGF=ECMandGFs)
repCV<-data.frame(rep=1:4, stringsAsFactors=FALSE)

#TODO: try to figure out the number of replicates and channels from the dimensions
nrReplicates<-4
nrChannels<-3

for(ch in 1:nrChannels){
  for(r in 1:nrReplicates){#Go through all the replicate plates in the screen
    for(ECMandGF in ECMandGFs){#Get the name of one of the ECM+GF combinations
      #Create column names for the CVs, CV counts and median intensity of the ECM+GF combos for each replicate plate
      meanColumnName<-paste0("MeanIntentTechRep",r,"Ch",ch)
      sdColumnName<-paste0("sdIntentTechRep",r,"Ch",ch)
      CVColumnName<-paste0("CVofTechRep",r,"Ch",ch)
      CVColumnCount<-paste0("CountofCVofRep",r,"Ch",ch)
      
      #Put the mean of the normalized values of the current channel,replicate and ECMandGF
      #into a new MeanTechRepXChX column
      ECMandGFDF[ECMandGFDF$ECMandGF==ECMandGF,meanColumnName]<-
        mean(ttReps[ttReps$ECMandGF==ECMandGF,paste0("normalized_r",r,"_ch",ch)],na.rm=TRUE)
      # Put the sd of the normalized values of the current channel and current ECM+GF
      # into a new sdTechRepXChX column
      ECMandGFDF[ECMandGFDF$ECMandGF==ECMandGF,sdColumnName]<-
        sd(ttReps[ttReps$ECMandGF==ECMandGF,paste0("normalized_r",r,"_ch",ch)],na.rm=TRUE)
      #Put the count of the non NA values of the current channel,replicate and ECMandGF
      #into a new CountofCVofRepXChX column
      ECMandGFDF[ECMandGFDF$ECMandGF==ECMandGF,CVColumnCount]<-
        sum(!is.na(ttReps[ttReps$ECMandGF==ECMandGF,paste0("normalized_r",r,"_ch",ch)]))
      
    }#end for name in tt$ECMandGF
    #Put the CV of the normalized values of the current channel and current ECM+GF
    #into a new CVofTechRepXChX column
    ECMandGFDF[,CVColumnName]<-ECMandGFDF[,sdColumnName]/abs(ECMandGFDF[,meanColumnName])
    
  }#End for r in 1:nrReplicates
}#End for ch in 1:nrChannels

for(ch in 1:nrChannels){
  for (ECMandGF in ECMandGFs){#Add the biological counts to the ECM and GF dataframe
    #Create a output name for the Biological replicate column
    intentColumnName<-paste0("BioRepCh",ch)
    #Create the names for the input columns that hold the technical replicate data
    intentDataColumnName<-paste0("MeanIntentTechRep",1:nrReplicates,"Ch",ch)
    #Calculate biological replicate values as the median intensity of the technical replicates
    #Put the biological replicate value into a new column for the current channel
    ECMandGFDF[ECMandGFDF$ECMandGF==ECMandGF,intentColumnName]<-
      median(as.numeric(ECMandGFDF[ECMandGFDF$ECMandGF==ECMandGF,intentDataColumnName]), na.rm=TRUE)
    #Create the output name for the biological replicate mean column
    meanColumnName<-paste0("BioRepMeanCh",ch)
    #Create the output name for the biological replicate sd column
    sdColumnName<-paste0("BioRepSDCh",ch)
    #Create the names for the input columns that hold the normalized intensities
    intentDataColumnName<-unlist(paste0("normalized_r",1:nrReplicates,"_ch",ch))
    #Calculate the mean of all replicates of the current ECMandGF
    ECMandGFDF[ECMandGFDF$ECMandGF==ECMandGF,meanColumnName]<-
      mean(unlist(ttReps[ttReps$ECMandGF==ECMandGF,intentDataColumnName]),na.rm=TRUE)
    #Calculate the sd of all replicates of the current ECMandGF
    ECMandGFDF[ECMandGFDF$ECMandGF==ECMandGF,sdColumnName]<-
      sd(unlist(ttReps[ttReps$ECMandGF==ECMandGF,intentDataColumnName]),na.rm=TRUE)
  }#End ECMandGFs for the current channel and replicate
  #Create the output name for the biological replicate CV column
  CVColumnName<-paste0("BioRepCVCh",ch)
  #Calculate the CV of all replicates of the current for the current channel
  ECMandGFDF[,CVColumnName]<-ECMandGFDF[,sdColumnName]/abs(ECMandGFDF[,meanColumnName])
}#End all replicates for ECMandGF CVs

for (ECMandGF in ECMandGFs){#Put the count data into a new column
  #Create the output name for the biological replicate count
  countColumnName<-paste0("CountofBioRepSpots")
  #Create the names for the input columns of channel 3 that hold the technical replicate count data
  countDataColumnName<-paste0("CountofCVofRep",1:nrReplicates,"Ch3")
  #Put the count of the non-NA values that went into the biological reps into a new column
  ECMandGFDF[ECMandGFDF$ECMandGF==ECMandGF,countColumnName]<-
    sum(ECMandGFDF[ECMandGFDF$ECMandGF==ECMandGF,countDataColumnName])
}#End of biological rep counts for all ECMandGFs 

#Calculate CVs for each channel of each replicate based on the techReps of ttReps
for(r in 1:nrReplicates){#create the dataframe of all CVs for the ECM and GF combinations
  #technical replicates of all replicate plates
  for(ch in 1:nrChannels){
    #calculate the per replicate per channel CVs 
    repName<-paste0("rep",r)
    columnName<-paste0("ch",ch)
    #Create an entry with the CVs
    dataColumnName<-paste0("techRep",r,"CV","Ch",ch)
    repCV[repCV$rep==r,columnName]<-median(ttReps[,which(names(ttReps)==paste0("techRep",r,"CV","Ch",ch))],
                                           na.rm=TRUE)
  }#End channels for the current replicate
}#End all replicates for ECMandGF CVs
#################End calculations and begin plots


#Open a pdf file to capture the plots
pdf(file=paste0(simpleName,".pdf"))
layout(matrix(c(1:12), 4, 3, byrow = TRUE))

#Plot histograms of the CVs for each channel of each replicate based of the ECM and GF combination
for(r in 1:nrReplicates){#create the dataframe of all CVs for the ECM and GF combinations
  #technical replicates of all replicate plates
  for(ch in 1:nrChannels){
    columnName<-paste0("CVofTechRep",r,"Ch",ch)
    CVs<-ECMandGFDF[ECMandGFDF[,columnName]<5,columnName]
    hist(x=CVs, main=paste0("ECM and GF CVs of Replicate ",r," Channel ",ch),
         breaks = seq(from=0,to=5,by=.1),ylim=c(0,200),cex.main=.8)
    abline(v=0.2,col="blue")
  }#End channels for the current replicate
}#End all replicates for ECMandGF CVs histograms

#Plot histograms of the filtered CVs for each channel of each replicate based of the ECM and GF combination
for(r in 1:nrReplicates){#create the dataframe of all CVs for the ECM and GF combinations
  #technical replicates of all replicate plates
  for(ch in 1:nrChannels){
    columnName<-paste0("CVofTechRep",r,"Ch",ch)
    countColumnName<-paste0("CountofCVofRep",r,"Ch",ch)
    CVs<-ECMandGFDF[ECMandGFDF[,countColumnName]>3,columnName]
    hist(x=CVs, main=paste0("Filtered ECM and GF CVs of Replicate ",r," Channel ",ch),
         breaks = seq(from=0,to=5,by=.1),ylim=c(0,200),cex.main=.8)
    abline(v=0.2,col="blue")
  }#End channels for the current replicate
}#End all replicates for ECMandGF CVs histograms

for(ch in 1:nrChannels){#Plot histograms of the biological replicate intensity values
  columnName<-paste0("BioRepCh",ch)
  normalizedIntensities<-ECMandGFDF[,columnName]
  hist(x=normalizedIntensities,
       main=paste0("ECM and GF Biological Replicate Intensities\nChannel ",ch),
       breaks = seq(from=0,to=8,by=.2),ylim=c(0,200),cex.main=.8)
}#End channels for the current replicate

for(ch in 1:nrChannels){#Plot histograms of the biological replicate CV values
  columnName<-paste0("BioRepCVCh",ch)
  CVs<-ECMandGFDF[,columnName]
  hist(x=CVs,
       main=paste0("ECM and GF Biological Replicate Intensities\nChannel ",ch),
       breaks = seq(from=0,to=5,by=.1),ylim=c(0,200),cex.main=.8)
}#End channels for the current replicate

#Plot the histogram of the count data
countColumnName<-"CountofBioRepSpots"
biologicalCounts<-ECMandGFDF[,countColumnName]
hist(x=biologicalCounts,
     main=paste0("ECM and GF Biological Replicate Counts"),
     breaks = seq(from=0,to=24,by=2),ylim=c(0,125),cex.main=.8)
abline(v=median(biologicalCounts), col="blue")

plot.new()

dev.off()#close the pdf file

#Write the ECMandGF dataframe to disk
write.table(x=ECMandGFDF,file=paste0("ECMandGF",simpleName,".txt"),quote=FALSE,sep="\t",row.names=FALSE)

