#Analyze MEA microarrays
#The input files from the Array-Pro are per microarray with one sheet per subarray
#There are four replicates of each plate treatment
#There are 384 combinations of 16 ECM proteins and 24 Growth Factors
#Each micrarray has 6 replicate treatements.
#The ECM/GF combinations are highly correlated with position on the microarray 

#Version 2.1
#Fixed flag for spots with few or no cells. It is treating replicate 1 different from the others
#low spots in plate one only are being "flagged". All low spots are changed to NAs

#Version 2.2
#Reduced analysis to Cellmask channel

#Version 2.3
#Changed name to MEA, putting functions into an environment and pushing to GitHub

#Author: Mark Dane


#Read the CMA functions into an environment
CMAEnv <- new.env()
sys.source("./Code/CMAFunctions.R", envir = CMAEnv)
attach(CMAEnv)

library(gdata) #read xls file
library(cellHTS2)#well plate pipeline functions
library(limma)#read GAL file

scriptName<-'./Code/MEA.R'
#Choose between Background, Raw or Net
chState<-"Net"
pinTipMethod<-'none'
#Do not use the log of the data values
logInts<-FALSE
#set the plate normalization method
normMethod<-'median'
#Don't have normalizePlates log the data
npLog<-FALSE
#Name the working directory and locate the data and metadata files
pathData<-paste("MEA Data")
#Create a list of the datafile names
dataFiles<-dir(path=pathData,full.names=TRUE,pattern=".xlsx")
pathGAL<-"~/Documents/Sophia MEMa 2/Data/"
GALfile<-list.files(path=pathData,pattern="gal")

cellTypes<-c("10A")
experimentName<-paste(cellTypes,", Input ",chState, ", pin ", pinTipMethod,", plate ",normMethod, sep="")

#Read the GAL file to get the data and layout
gal<-readGAL(GALfile,path=pathData)
layout<-getLayout(gal)
#Delete the Annotation column
gal<-subset(gal,select=-c(Annotation))
#Create a plate column
gal$plate<-1
#Rename the GAL columns
names(gal)<-c("Grid","Row","Column","ECM","GrowthFactor","plate")
#Create a column that is the combination of the ECM and GF names
gal$name<-paste(gal$ECM,gal$GrowthFactor)
#Add arrayRow, arrayColumn and PlateIndex columns and reorganize gal to read 
#across rows then down columns
gal<-orderDFV2(gal,gridsPerRow=layout$ngrid.c)
#Read each datafile and merge it with the gal annotations
dfs<-lapply(X=dataFiles,FUN=readReplicates,gal)
#Simplify the channel names
dfs<-lapply(X=dfs,FUN=simplifyNames)
#Make the data values numeric
dfs<-lapply(X=dfs,FUN=makeNumeric)
#Rbind the dataframes from each replicate into one dataframe 
inputDF<-rbindDfs(dfs)
#Remove the 488 and 532 channel data. This is done to simplify the dataset.
inputDF<-subset(inputDF,select=-grep(pattern="488|532",names(inputDF)))
#Assign well alhpanumeric names
inputDF$Well<-with(inputDF,convertWellCoordinates(x=PlateIndex,
                              pdim=c("nrow"=max(arrayRow),"ncol"=max(arrayColumn)))$letnum)
#Add a controls type column
inputDF<-addControls(inputDF,pos=c("Col I PBS"),neg=c("MG PBS"))
#Build the raw object using the chState to select the correct Teacan column types
#The nested names and grep code causes the channels to use the name from inputDF
x <- buildCellHTS2(data.frame(inputDF[names(inputDF[grep(pattern=chState,names(inputDF))])],
                              plate=inputDF$plate,replicate=inputDF$replicate, well=inputDF$Well,
                              stringsAsFactors=FALSE))
name(x)<-experimentName
par(mfrow=c(3,1))

#Filter out low cell spots
x<-filterOnCellCount(df=x,countCh=1,thresh=50)
#Log the values if the logInts flag was set
if(logInts){
  x<-logObject(x, base=2)
}

#Configure the object.
x<-configure(x,confFile=(function() list(as.integer(c(1,dim(x)), length=2),
                data.frame(cbind(Plate=inputDF$plate, Well=inputDF$Well,
                                 Content=inputDF$controlStatus),
                           stringsAsFactors=FALSE))),  path=pathData, descripFile='/MEA Description.txt')
#Annotate the raw object
x<-chtsAnnotate(annDF=data.frame(Plate=inputDF$plate, Well=inputDF$Well, ECM=inputDF$ECM,
                                 GrowthFactor=inputDF$GrowthFactor,
                                 GeneID=inputDF$name,Grid=inputDF$Grid, Row=inputDF$Row,
                                 Column=inputDF$Column)[1:(dim(inputDF)[1]/dim(x)[2]),],x=x)
#Get the spatial variation indices for the raw data
SVI<-spatialVariation(x=normalizePlates(x,scale="multplicative",
                                        log=FALSE,
                                        method='locfit',save.model=TRUE,
                                        varianceAdjust="none", nn=.01))

#Write the raw data report
out <- writeReport(x, outdir=paste(experimentName, "raw"),force=TRUE,
                   mainScriptFile=scriptName,
                   settings=list(plateList=list(reproducibility=list(include=FALSE, map=FALSE),
                                intensities=list(include=TRUE, map=FALSE,thumbFactor=.5,
                                range=quantile(Data(x), probs=c(.01,.95), na.rm=TRUE))),
                                plateConfiguration=list(include=FALSE)))
#Normalize the intensities 
xn<-chtsNormalizePlates(x,normMethod='median',logged=logInts)
xsc<-summarizeReplicates(xn, summary="median")
#Write the HTML report of the normalized data
out<-writeReport(raw=x, normalized=xn, scored=xsc,
                 force = TRUE, outdir = paste(experimentName),mainScriptFile=scriptName,
                 settings=list(plateList=list(reproducibility=list(include=FALSE, map=FALSE),
                          intensities=list(include=TRUE, map=FALSE,thumbFactor=.5,
                          range=quantile(Data(xn), probs=c(.01,.95),na.rm=TRUE))),
                          plateConfiguration=list(include=FALSE),
                          screenSummary=list(scores=list(range=c(-3,3), map=FALSE))))

#Plot the pingrid boxplots of the raw and normalized data
par(mfrow=c(2,2),mar=c(4,4,2,2))
chtsPinBoxPlots(x=x,main="",col="pink",ylim=c(50,40000))
chtsPinBoxPlots(x=xn,main=" ",col="light blue")

detach(name=CMAEnv)
