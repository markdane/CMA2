#Create and analyze the Simulated dataset
#Author: Mark Dane


#Version 1.7

library(gdata)    #read xls file
library(cellHTS2) #well plate pipeline functions
library(limma)    #read GAL file
library(pROC)     #for ROC analysis

#Read the CMA functions into an environment
CMAEnv <- new.env()
sys.source("./Code/CMAFunctions.R", envir = CMAEnv)
sys.source("./Code/SimulateFunctions.R", envir=CMAEnv)
attach(CMAEnv)

scriptName<-'./Code/Simulated.R'

#Select the parameter and method packet that sets global variables for 
#creating the simulated dataset
source('./Code/plateEllLine.R')

#Build the dataframe of simulated data 
base<-simulatePlates(nrPlates=nrPlates,nrChannels=nrChannels,
                     plateGridRows=plateGridRows,plateGridCols=plateGridCols,gridRows=gridRows,gridCols=gridCols,
                     plateMeans=plateMeans,plateSDs=plateSDs,
                     pinMeans=pinMeans,pinSDs=pinSDs,
                     channelMeans=channelMeans,channelSDs=channelSDs,
                     pctGenesOfInterest=pctGenesOfInterest, ellipList=ellipList, lineList=lineList)

#Use the base plate to make the replicate plates
xd<-makeReplicates(base, replicateMeans=replicateMeans,replicateSDs=replicateSDs)

#Add spatial perturbations to the choosen channel
xd$channel1[xd$replicate==1]<-spatPertPlate(plate=xd$channel1[xd$replicate==1],plateRows=plateRows,plateCols=plateCols,
                                            ellipList=ellipList,lineList=lineList)
#Choose the bivariate local regression normalization method
normMethod="locfit"
nn=.01

#Set the output parameters for Write Report
setSettings(list(plateList=list(reproducibility=list(include=FALSE, map=TRUE), intensities=list(include=TRUE, map=TRUE, thumbFactor=.3)),
                 screenSummary=list(scores=list(range=c(-3,3), map=TRUE))))

#Build the raw object based on the xd dataframe
x <- buildCellHTS2(xd[grep(pattern="channel|well|plate|replicate",x=names(xd))])
name(x)<-experimentName

#Filter out low cell spots
x<-filterOnCellCount(df=x,countCh=1,thresh=0)

if(logInts){
  x<-logObject(x, base=2)
  plateNormScale<-'additive'
}

#Configure the object as all samples
x<-configure(x, confFile=(function() list(as.integer(c(nrPlates,plateRows*plateCols), length=2),
                                    data.frame(cbind(Plate=rep(xd$plate,times=nrReplicates),
                                    Well=xd$well, Content=xd$content),stringsAsFactors=FALSE))),
             path="Simulated Data", descripFile=descripFile)
#Create annotations that include target names and the block, row and column position
x<-chtsAnnotate(annDF=data.frame(Plate=xd$plate, Well=xd$well,GeneID=paste('Gene',xd$well, sep=""),
                      Grid=xd$grid,Row=xd$row,Column=xd$column)[1:(dim(xd)[1]/dim(x)[2]),],x=x)

out <- writeReport(x, outdir=paste("Raw", experimentName),force=TRUE, mainScriptFile=scriptName,
                   settings=list(plateList=list(reproducibility=list(include=FALSE, map=TRUE),
                            intensities=list(include=TRUE, map=TRUE,thumbFactor=.3)),
                                 screenSummary=list(scores=list(range=c(-3,3), map=TRUE))))

#Locfit normalize at the microarray level
xpln<-chtsNormalizePlates(x,normMethod='locfit',logged=logInts,nn=nn)

xsc<-summarizeReplicates(xpln, summary="median")

#Get the spatial variation indices for the raw data
SVI<-spatialVariation(x=normalizePlates(x,scale="additive",
                                        log=FALSE,
                                        method='locfit',save.model=TRUE,
                                        varianceAdjust="none", nn=nn))


#Write the HTML report of the normalized data
out<-writeReport(raw=x, normalized=xpln, scored=xsc,
                 force = TRUE, outdir = paste("Normalized",experimentName),mainScriptFile=scriptName,
                 settings=list(plateList=list(reproducibility=list(include=FALSE, map=FALSE),
                          intensities=list(include=TRUE, map=FALSE,thumbFactor=.5,
                          range=quantile(Data(xpln), probs=c(.01,.95),na.rm=TRUE))),
                              plateConfiguration=list(include=FALSE),
                               screenSummary=list(scores=list(range=c(-3,3), map=FALSE))))

########## Median normalizing and scoring for baseline ########
#Create a median normalized object as a baseline
xMedian<-chtsNormalizePlates(x,normMethod='median',logged=logInts)

#Summarize and score the median normalized object
xMediansc<-scoreReplicates(xMedian, sign="+", method="zscore")
xMediansc<-summarizeReplicates(xMediansc, summary="mean")
#Get a topTable from the median normalized object
outMedian<-getTopTable(cellHTSlist=list("raw"=x, "normalized"=xMedian,
                                        "scored"=xMediansc), file=tempfile())
#Get the ppvs after median normalization
rocobjMedian<-roc(response=(outMedian$wellAnno=="other"|outMedian$wellAnno=="pos"),
                  predictor=.5*(outMedian$normalized_r1_ch1+outMedian$normalized_r2_ch1))
mppv<-coords(rocobjMedian ,ret="ppv",x="best")
#Get the ppvs using the median of the raw replicate values
rocobjRaw<-roc(response=(outMedian$wellAnno=="other"|outMedian$wellAnno=="pos"),
            predictor=(outMedian$median_ch1))
rppv<-coords(rocobjRaw ,ret="ppv",x="best")
#Start the ROC plot with the raw and median curves
par(mfrow=c(1,1))
plot(rocobjRaw,xlim=c(1,.8),ylim=c(.7,1),
     col="pink",lty=2, lwd=5)
plot(rocobjMedian,add=TRUE,col="light blue",lty=3,lwd=5)
legend(x="bottomright",legend=c("raw","median normalized","spatially normalized \n nn = [.5, .001]"),
       fill=c("pink","light blue","purple"))

########## End Median normalizing and scoring for baseline ########

nnS<-c(.5,.1,.05,.01,.005,.001)
nppvS<-NULL
zppvS<-NULL
RPppvS<-NULL

#Create a series of ROC results for different locfit nn values
for (nn in nnS){
  normMethod="locfit"
  pinTipMethod="none"
  
  #Locfit normalize at the microarray level
  xpln<-chtsNormalizePlates(x,normMethod='locfit',logged=logInts,nn=nn)
  
  #Get the ppv based on mean spatial normalization values
  rocobjN<-roc(response=(wellAnno(xpln)=="other"|wellAnno(xpln)=="pos"),
               predictor=apply(Data(xpln),MARGIN=1,FUN=sum)/dim(Data(xpln))[2])
  ppv<-coords(rocobjN ,ret="ppv",x="best")
  nppvS<-c(nppvS, ppv)
  plot(rocobjN,add=TRUE,col="purple",lty=2)
  
  #Perform z score scoring
  xsc<-scoreReplicates(xpln, sign="+", method="zscore")
  xsc<-summarizeReplicates(xsc, summary="mean")
  
  #Perform rank product scoring
  xRP<-chtsRankReplicates(x=xpln,channels="channel1")
  
  #Get the ppv after locfit normalization based on rank product score
  rocobjRP<-roc(response=(wellAnno(xRP)=="other"|wellAnno(xRP)=="pos"),
              predictor=Data(xRP))
  ppv<-coords(rocobjRP ,ret="ppv",x="best")
  RPppvS<-c(RPppvS, ppv)
  
  #Get the ppv after normalization based on z score
  rocobjZ<-roc(response=(wellAnno(xsc)=="other"|wellAnno(xsc)=="pos"),
              predictor=Data(xsc))
  ppv<-coords(rocobjZ ,ret="ppv",x="best")
  zppvS<-c(zppvS, ppv)
}#End loop of nn values

#Create Rank Product vs Z score graph at all nn values
plot(y=zppvS,x=log10(nnS),col="blue", pch=19, type="l",ylim=c(.3,1),lwd=3,
     #main="Positive Predictive Rate vs Normalization Parameter nn",
     ylab="Positive Predictive Value",xlab="Log10 of nn")
points(x=log10(nnS),y=RPppvS, type="l",lwd=3,col="violetred")
legend(x="topright",legend=c("z scored","rank product scores"),
       fill=c("blue","violetred"))

plot(y=nppvS,x=log10(nnS),col="purple", pch=19, type="l",ylim=c(.3,1),lwd=3,axes=TRUE,
     # main="Positive Predictive Rate vs Normalization Parameter nn",
     ylab="Positive Predictive Rate",xlab="Log10 of nn")
abline(h=rppv, col="pink",lwd=3)
abline(h=mppv, col="light blue",lwd=3,lty=2)
legend(x="topright",legend=c("raw","median normalized","spatially normalized"),
       fill=c("pink","light blue","purple"))

########Evaluate PPV for Z score vs Rank Product at different list sizes
plot(x=NULL,xlim=c(0,350),ylim=c(0,1),ylab="Positive Predictive Value",xlab="List Size")
legend(x="topright",legend=c("z scored (nn)","rank product scored"),
       fill=c("blue","purple"))
nnS<-c(1,.5,.1)
for (nn in nnS){
#Locfit normalize at the microarray level
xpln<-chtsNormalizePlates(x,normMethod='locfit',logged=logInts,nn=nn)
#Perform z score scoring and order by the Z scores
xsc<-scoreReplicates(xpln, sign="+", method="zscore")
xsc<-summarizeReplicates(xsc, summary="mean")
xscD<-data.frame(score=Data(xsc),wellAnno=wellAnno(xsc), stringsAsFactors=FALSE)
xscD<-xscD[order(xscD[1]),]
#Perform rank product scoring and order by the rank products
xRP<-chtsRankReplicates(x=xpln,channels="channel1")
xRPD<-data.frame(score=Data(xRP),wellAnno=wellAnno(xRP), stringsAsFactors=FALSE)
xRPD<-xRPD[order(xRPD[1]),]

ZppvList<-NULL
RPppvList<-NULL
for(s in 1:350){
      #compare z score to rank product at various list sizes
      ZppvList<-c(ZppvList,sum(xscD[1:s,'wellAnno']%in% c("pos","other","flagged"))/s)
      RPppvList<-c(RPppvList,sum(xRPD[1:s,'wellAnno']%in% c("pos","other","flagged"))/s)
}
points(ZppvList,type='l',col='blue',lwd=3)
text(x=100,y=ZppvList[100]-.035,labels=nn)
points(RPppvList,type='l',col='purple',lwd=3)
}#End figure for z scored and rank product at different nn and list sizes

######### rank product scoring ###########

# #Create Rank Product vs Z score graphs at all nn values
# points(y=zppvS,x=log10(nnS),col="blue", pch=19, type="l",ylim=c(.3,1),lwd=3,
#      #main="Positive Predictive Rate vs Normalization Parameter nn",
#      ylab="Positive Predictive Value",xlab="Log10 of nn")
# points(x=log10(nnS),y=RPppvS, type="l",lwd=3,col="violetred")
# legend(x="topright",legend=c("z scored","rank product scores"),
#        fill=c("blue","violetred"))
# coords.roc(rocobj, "best", ret=c("sens", "spec", "ppv", "npv","tp","fp"))
# 
# ######### end rank product scoring ###########


detach(name=CMAEnv)
