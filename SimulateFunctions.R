#Functions that create simulated data for cell-based microarrays
#Author: Mark Dane

#Version 1.7
#Reformatted to a fixed name 


fNormIntensities<-function(mean,sd,nrRows,nrColumns){
  #Return a vector of folded normally distributed intensities
  #Fold back any negative values to simulate intensities > 0
  #Get intensities for the grid based on grid mean and variance, then modify them by the pin mean
  set.seed(1234)
  return(abs(rnorm(nrRows*nrColumns, mean=mean, sd=sd)))
}#End fNormIntensities

assignGrids<-function(gridsPerRow,gridsPerCol,gridRows,gridCols){
  #Create the array of pin grid,row and column assignments
  #####Grid assignments #######
  #Create a sequence of the first pin grid in each pin grid row
  firstGrid<-seq(from=1,to=gridsPerCol*gridsPerRow,by=gridsPerRow)
  
  #create a list of matrices with the pin grid numbers of each row
  gridRow<-lapply(firstGrid,FUN=function(x,gridCols,gridRows,gridsPerRow){
    matrix(data=rep(x:(x+gridsPerRow-1),each=gridCols,times=gridRows),
           nrow=gridRows,byrow=TRUE)},
    gridCols=gridCols,gridRows=gridRows,gridsPerRow=gridsPerRow)
  #Make one matrix out of all of the grid numbers
  grids<-do.call("rbind", gridRow)
  
  ######Row assignments ############
  #Create the matrix of pin grid row assignments
  rows<-matrix(data=rep(x=1:gridRows,each=gridCols),nrow=gridRows,byrow=TRUE)
  
  #Create one complete microarray row of the pin rows
  gridRowsRow<-integer()
  for(i in 1:gridsPerRow) gridRowsRow<-cbind(gridRowsRow,rows)
  #Create the complete microarray pin row assignments
  maGridRows<-integer()
  for(i in 1:gridsPerCol) maGridRows<-rbind(maGridRows,gridRowsRow)
  
  #######Column Assignments #########
  #Create the matrix of pin grid column assignments
  cols<-matrix(data=rep(x=1:gridCols,times=gridRows),nrow=gridRows,byrow=TRUE)
  
  #Create one complete microarray row of the pin rows
  gridColsRow<-integer()
  for(i in 1:gridsPerRow) gridColsRow<-cbind(gridColsRow,cols)
  #Create the complete microarray pin row assignments
  maGridCols<-integer()
  for(i in 1:gridsPerCol) maGridCols<-rbind(maGridCols,gridColsRow)
  #########
  
  #return an array that reads out across the rows of the microarray
  return(list(array(t(grids)),array(t(maGridRows)),array(t(maGridCols))))
  
}#End assignGrids function


gridIntensities<-function(plateGridRows,plateGridCols,gridRows,gridCols,plateMean,plateSD,
                          pinMeans,pinSDs,channelMean,channelSD){
  #Create one channel of intensities, one plate, one grid at a time
  plateRows<-plateGridRows*gridRows
  plateCols<-plateGridCols*gridCols
  #Create an empty matrix for the intensities
  plateIntensities<-matrix(nrow=plateRows, ncol=plateCols)
  index<-1
  #step down the rows of the grids in a plate
  for(plateGridRow in 1:plateGridRows){
    #step across the columns of grids in a plate
    for(plateGridCol in 1:plateGridCols){
      #Find the 1,1 corner of the next grid in plate coordinates
      plateCorner<-c((plateGridRow-1)*gridRows+1,(plateGridCol-1)*gridCols+1)
      #use the corner to place the grid values
      #Create normally distributed values, this can be modified to give each
      #grid a different distribution, mean, etc.
      mean<-plateMean*pinMeans[index]*channelMean
      sd<-plateSD+pinSDs[index]+channelSD
      plateIntensities[plateCorner[1]:(plateCorner[1]+(gridRows-1)),
                       plateCorner[2]:(plateCorner[2]+(gridCols-1))]<-
        fNormIntensities(mean=mean,sd=sd,nrRows=gridRows,nrColumns=gridCols)
      #Keep track of the blocks being created by creating a matrix of block indices
      index<-index+1
    }
  }
  #convert from matrix to array so that it can become a column in the dataframe
  plateIntensities<-as.numeric(t(plateIntensities)) 
  return(plateIntensities)
}

configArray<-function(xd,negPositions=c(3,8,8,3),posPositions=c(3,3,8,8,4,6)){
  #Configure the dataframe as samples, positive and negative controls.
  #Use the grid, row and block indices. All grids, plates and replicates will have 
  #the same controls configuration
  #TODO: Generalize this function for the number of each type of control
  
  #Start with all wells configured as samples
  xd$content<-'sample'
  #Put positive controls in the positions passed through. Use the first element
  #as the row and the next as the column. 
  xd$content[xd$row==posPositions[1]&xd$column==posPositions[2]]<-'pos'
  xd$content[xd$row==posPositions[3]&xd$column==posPositions[4]]<-'pos'
  xd$content[xd$row==posPositions[5]&xd$column==posPositions[6]]<-'pos'
  xd$content[xd$row==negPositions[1]&xd$column==negPositions[2]]<-'neg'
  xd$content[xd$row==negPositions[3]&xd$column==negPositions[4]]<-'neg'
  return(xd)
}

ePerturb<-function(e, plate, plateRows, plateCols){
  #Perturb the plate intensities using the described ellipsoidal
  #use the paramerizations of an ellipsoid to define the perturbations
  #xc is the x (column) location of the center
  #yc is the y (row) location of the center
  #ax is the radius along the x (column) axis
  #ay is the radius along the y (row) axis
  #p is the highest amount of perturbation
  xc<-e[1]
  yc<-e[2]
  ax<-e[3]
  ay<-e[4]
  p<-e[5]
  #create a matrix that contains the perturbation factors
  pMat<-matrix(nrow=plateRows, ncol=plateCols)
  for (row in 1:plateRows){
    for (col in 1:plateCols){
      #Calculate perturbations for the ellipse
      zsquared<-p^2*(1-(col-xc)^2/ax^2-(row-yc)^2/ay^2)
      if(zsquared>=0){
        z<-sqrt(zsquared)}
      else z<-0
      if (p>0){
      pMat[row,col]<-1+z}
      else pMat[row,col]<-1-z
    }#End column perturbations
  }#End row perturbations
  
  #multiply the input plate intensities by the perturbations
  return(plate<-plate*as.numeric(t(pMat)))
}#End ePerturb function


lPerturb<-function(pts, plate, plateRows, plateCols){
  #Perturb the plate intensities within a rectangle between the given points
  
  #x1,y1  x2,y2 endpoints of the line
  #t thickness in spots of the line
  #p amount of perturbation, added to 1 and then multiplied
  x1<-pts[1]
  y1<-pts[2]
  x2<-pts[3]
  y2<-pts[4]
  t<-pts[5]
  p<-pts[6]

  #create a matrix that contains the perturbation factors
  pMat<-matrix(nrow=plateRows, ncol=plateCols)
  for (row in 1:plateRows){
    for (col in 1:plateCols){
      #use the equation of the line to check if within t distance from it
      if((x2-x1)!=0) m<-(y2-y1)/(x2-x1)
      else m=999
      b<-y1-(m*x1)
      d<-abs(m*row-col+b)/sqrt(m^2+1)
      if ((d<=t)&(row>=min(x1,x2))&(row<=max(x1,x2))&(col>=min(y1,y2))&(col<=max(y1,y2))) pMat[row,col]<-1+p
      else pMat[row,col]<-1
    }#End column perturbations
  }#End row perturbations
  #multiply the input plate intensities by the perturbations
  return(plate<-plate*as.numeric(t(pMat)))
}#End lPerturb function

simulatePlates<-function(nrPlates, nrChannels, plateGridRows, plateGridCols, gridRows, gridCols,
                         plateMeans, plateSDs, pinMeans,pinSDs,
                         channelMeans,channelSDs,pctGenesOfInterest, ellipList, lineList){
  #Create and return a dataframe of intensity values based on the input parameters
  
  nrPinGrids<-plateGridRows*plateGridCols
  
  #Calculate the number of rows, columns and wells in a plate based on the grid information
  plateRows<-plateGridRows*gridRows
  plateCols<-plateGridCols*gridCols
  nrWells<-plateRows*plateCols
  
  #Create well names to match the size of the plate
  wells<-convertWellCoordinates(x=1:nrWells,pdim=c("nrow"=plateRows,"ncol"=plateCols))$letnum
  #Create a dataframe for all replicates, plates and channels
  #The rows in the dataframe will be the number of rows in one plate x number of plates x number of replicates
  xd <- expand.grid(well=wells, plate=1:nrPlates, stringsAsFactors = FALSE)
  rm(wells)
  
  #Get the grid indices (grid, row and column) in plate coordinates and bind as columns in the dataframe
  gridList<-assignGrids(gridsPerRow=plateGridCols,gridsPerCol=plateGridRows,
                        gridRows=gridRows,gridCols=gridCols)
  
  xd<-cbind(xd,grid=gridList[[1]],row=gridList[[2]],column=gridList[[3]])
  
  #Create a character vector for the content configuration
  xd<-configArray(xd,negPositions=c(3,8,8,3),posPositions=c(3,3,8,8,4,6))
  
  #Get all channels of all plates in the first and maybe only replicate set
  intensities<-numeric(length=dim(xd)[1])
  for(ch in 1:nrChannels){
    for(p in 1:nrPlates){
      start<-(p-1)*nrWells+1
      intensities[start:(start+nrWells-1)]<-
        abs(gridIntensities(plateGridRows=plateGridRows,plateGridCols=plateGridCols,
                        gridRows=gridRows,gridCols=gridCols,plateMean=plateMeans[p],
                        plateSD=plateSDs[p], pinMeans=pinMeans,pinSDs=pinSDs,
                        channelMean=channelMeans[ch],channelSD=channelSDs[ch]))
    }
    #Add the first replicate set of plate values to the dataframe and name it
    xd<-cbind(xd,intensities)
    names(xd)[dim(xd)[2]]<-paste("channel",ch, sep="")
  }
  #Select a % of the samples to be the genes of interest
  genesOfInterest<-xd$well[sample((1:nrWells)[xd$content[1:nrWells]=='sample'], nrWells*pctGenesOfInterest, replace=FALSE)]
  
  #Knockdown the genes of interest by a factor of geneMean
  set.seed(1234)
  GOIFactors<-abs(rnorm(length(genesOfInterest),mean=(1-GOIFactor), sd=GOISD))
  #knockdown all of the channels using the genes of interest fatcors
  xd[xd$well %in% genesOfInterest,(dim(xd)[2]-nrChannels+1):(dim(xd)[2])]<-
    xd[xd$well %in% genesOfInterest,(dim(xd)[2]-nrChannels+1):(dim(xd)[2])]*GOIFactors
  
  #Change the content of xd to make the genes of interest as type 'other'
  xd[xd$well%in%genesOfInterest,'content']<-'other'
  
  #Add a column to xd showing the GOI knockdown factor. Use 1 for no factor (not a gene of interest)
  xd$GOIfactor<-1
  xd[xd$well%in%genesOfInterest,'GOIfactor']<-GOIFactors
  
  #Reduce the positive controls intensities
  set.seed(1234)
  xd[xd$content=='pos',(dim(xd)[2]-nrChannels):(dim(xd)[2]-1)]<-
    xd[xd$content=='pos',(dim(xd)[2]-nrChannels):(dim(xd)[2]-1)]*
    abs(rnorm(sum(xd$content=='pos'),mean=(1-posFactor), sd=posSD))  
  return(xd)
}

makeReplicates<-function(base,replicateMeans,replicateSDs){
  #Make replicates from the base dataframe with one subarray having one or more channels
  # using the given means and standard deviations

  nrReplicates<-length(replicateMeans)
  #There are eight non-channel columns in the base dataframe. Use this to set nrChannels
  nrChannels<-dim(base)[2]-7
  #Find out how many plates are in the base dataframe
  nrPlates<-max(base$plate)
  #Calculate the number of rows and columns in the grids and one plate
  gridRows<-max(base$row)
  gridCols<-max(base$column)
  nrBlocks<-max(base$grid)
  nrWells<-nrBlocks*gridRows*gridCols
  
  #Make a copy of base in order to make the replicates
  repBase<-base
  for(r in 1:nrReplicates){  
    set.seed(1234)
    #Create the new channel values
    perturbs<-data.frame(matrix(abs(rnorm(nrWells*nrPlates*nrChannels,replicateMeans[r], replicateSDs[r])),
                                nrow=dim(base)[1]))
    repChannels<-perturbs*base[7:(6+nrChannels)]
    repBase<-data.frame(base[1:2],replicate=r,base[3:6],repChannels,base[dim(base)[2]])
    names(repBase)<-c(names(base)[1:2],'replicate',names(base)[-1:-2])
    if(r==1) outputBase<-repBase
    else outputBase<-rbind(outputBase,repBase)
  }
  #detach(name)
  return(outputBase)
}
makePVIReplicates<-function(base,replicateMeans,replicateSDs){
  #Make replicates from the base dataframe with one subarray having one or more channels
  # using the given means and standard deviations
  
  nrReplicates<-length(replicateMeans)
  #There are eight non-channel columns in the base dataframe. Use this to set nrChannels
  nrChannels<-dim(base)[2]-7
  #Find out how many plates are in the base dataframe
  nrPlates<-max(base$plate)
  #Calculate the number of rows and columns in the grids and one plate
  gridRows<-max(base$row)
  gridCols<-max(base$column)
  nrBlocks<-max(base$grid)
  nrWells<-nrBlocks*gridRows*gridCols
  
  #Make a copy of base in order to make the replicates
  repBase<-base

  for(r in 1:nrReplicates){  
    set.seed(1234)
    #Create the new channel values
    perturbs<-data.frame(matrix(abs(rnorm(nrWells*nrPlates*nrChannels,replicateMeans[r], replicateSDs[r])),
                                nrow=dim(base)[1]))
    repChannels<-perturbs*base[7:(6+nrChannels)]
    repBase<-data.frame(base[1:2],replicate=r,base[3:6],repChannels,base[dim(base)[2]])
    names(repBase)<-c(names(base)[1:2],'replicate',names(base)[-1:-2])
    if(r==1) outputBase<-repBase
    else outputBase<-rbind(outputBase,repBase)
  }
  #detach(name)
  return(outputBase)
}#End makePVIReplicates

createConfig<-function(plateGridRows,gridRows,plateGridCols,gridCols){
  #Create a character vector for the content configuration
  configuration<-matrix('sample',nrow=plateGridRows*gridRows, ncol=plateGridCols*gridCols)
  
  #Put the positive and negative intensity controls in each grid
  for(plateGridRow in 1:plateGridRows){
    #step down the grids in a plate
    for(plateGridCol in 1:plateGridCols){
      #Find the 1,1 corner of the next grid
      plateCorner<-c((plateGridRow-1)*gridRows+1,(plateGridCol-1)*gridCols+1)
      #use the corner to place the positive controls in the 3,3 and 8,8 positions
      configuration[plateCorner[1]+2,plateCorner[2]+2]<-'pos'
      configuration[plateCorner[1]+7,plateCorner[2]+7]<-'pos'
      configuration[plateCorner[1]+2,plateCorner[2]+7]<-'neg'
      configuration[plateCorner[1]+7,plateCorner[2]+2]<-'neg'
    }
  }
  configArray<-as.character(t(configuration))
  return(configArray)
}

######## Spatial perturbations ################
spatPertPlate<-function(plate,plateRows,plateCols,ellipList=NULL,lineList=NULL){
  #Apply spatial perturbations to a plate of intensities
  #add ellipsoid and linear perturbations to the plate intensities
  for (e in ellipList) plate<-ePerturb(unlist(e),plate,plateRows,plateCols)
  for (l in lineList) plate<-lPerturb(unlist(l),plate,plateRows,plateCols)
  return(plate)
}
