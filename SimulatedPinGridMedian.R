#Create and analyze the Simulated dataset using pin grid median normalization
#Author: Mark Dane

#Version 1.0

library(gdata)    #read xls file
library(cellHTS2) #well plate pipeline functions
library(limma)    #read GAL file
library(pROC)     #for ROC analysis

#Read the CMA functions into an environment
CMAEnv <- new.env()
sys.source("./Code/CMAFunctions.R", envir = CMAEnv)
sys.source("./Code/SimulateFunctions.R", envir=CMAEnv)
attach(CMAEnv)

scriptName<-'./Code/SimulatedPinGridMedian.R'

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
#Do not normalize at the microarray level
normMethod="none"
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

#Normalize at the pin grid level
xptn<-pinTipNormalize(x)
xpln<-chtsNormalizePlates(xptn,normMethod='none',logged=logInts,nn=nn)

xsc<-summarizeReplicates(xpln, summary="median")

#Get the spatial variation indices for the raw data
SVI<-spatialVariation(x=normalizePlates(x,scale="additive",
                                        log=FALSE,
                                        method='locfit',save.model=TRUE,
                                        varianceAdjust="none", nn=nn))


#Write the HTML report of the normalized data
out<-writeReport(raw=x, normalized=xpln, scored=xsc,
                 force = TRUE, outdir = paste("Pin Grid Normalized",experimentName),mainScriptFile=scriptName,
                 settings=list(plateList=list(reproducibility=list(include=FALSE, map=FALSE),
                          intensities=list(include=TRUE, map=FALSE,thumbFactor=.5,
                          range=quantile(Data(xpln), probs=c(.01,.95),na.rm=TRUE))),
                              plateConfiguration=list(include=FALSE),
                               screenSummary=list(scores=list(range=c(-3,3), map=FALSE))))


detach(name=CMAEnv)
