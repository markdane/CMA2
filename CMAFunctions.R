#Pin block normalization functions
#Mark Dane
#Version 1

#Version 1.1 adding replicate handling


medianBlock<-function(x,drawPlots=FALSE){
  #A low level normalization function 
  #x is a matrix of intensity values
  #x is returned centered about its median
  x<-x-median(x, na.rm=TRUE)
  if(drawPlots){
    plot(x, main="Median Fit")
  }
  return(x)
}

loessBlock<-function(ints, drawPlots=FALSE, span=.1){
  #A low level normalization function
  #ints is a matrix of intensity values
  #order is a print order that may be used
  #the default is to order the values by row
  #The residuals of the actual values from the loess fitted values are returned.

  z<-as.numeric(t(ints))
  cols<-1:dim(ints)[2]
  rows<-1:dim(ints)[1]
  x<-rep(cols,times=(dim(ints)[1]))
  y<-rep(rows,each=(dim(ints)[2]))
  loess.model<-loess(z~x*y,na.action='na.exclude', span=span)
  residuals<-z-predict(loess.model)
   if(drawPlots){
     image(y=cols,x=rows,matrix(predict(loess.model),nrow=length(rows),byrow=TRUE))
#     plot(x[order(loessOrder)], main="Loess Fit")
#     lines(predict(loess.model,order)[order(loessOrder)], col='red')
#     abline(median(x, na.rm=TRUE),b=0, col="blue")
   }
  return(residuals)
}

loessPlate<-function(ints, drawPlots=FALSE, span=.1){
  #A low level normalization function
  #ints is a matrix of intensity values
  #The residuals of the actual values from the loess fitted values are returned.
  
  z<-as.numeric(t(ints))
  cols<-1:dim(ints)[2]
  rows<-1:dim(ints)[1]
  x<-rep(cols,times=(dim(ints)[1]))
  y<-rep(rows,each=(dim(ints)[2]))
  loess.model<-loess(z~x*y,na.action='na.exclude', span=span)
  residuals<-z-predict(loess.model)
  #   if(drawPlots){
  #     plot(x[order(loessOrder)], main="Loess Fit")
  #     lines(predict(loess.model,order)[order(loessOrder)], col='red')
  #     abline(median(x, na.rm=TRUE),b=0, col="blue")
  #   }
  return(residuals)
}


pinBlock<-function(x,order=1:length(x), method='median', drawPlots=FALSE){
  #A pin block normalization wrapper
  
  #x is a matrix to be normalized. x will be read in by row first.
  #order is vector that orders the matrix for some normalizations. For instance,
  # the print order may be used for the loess normalization.
  #y must be a numeric vector with the same length as x
  #method is a character specifying the normalization method. Allowed values are "median"
  # "loess"
  switch(EXPR = method, median=medianBlock(x,drawPlots), loess=loessBlock(x,order,drawPlots))
    }

refChannelNormalize<-function(x,refChannel, logged=FALSE){
  #Return an object with channels normalized by dividing by the reference channel
  #Read the number of channels from the object
  nrChannels<-dim(Data(x))[3] 
  for(ch in 1:nrChannels){
    if(ch!=refChannel){
      if(!logged) Data(x)[,,ch]<-Data(x)[,,ch]/(Data(x)[,,refChannel]+0.001)
      else Data(x)[,,ch]<-Data(x)[,,ch]-(Data(x)[,,refChannel])
    }
  }
  x@state["normalized"]<-TRUE
  return(x)
}

pinTipNormalize<-function(x,order=1:pdim(x)[1],method='median',drawPlots=FALSE){
  #pin tip normalize the cellHTS object x using method (only median and loess are supported)
  #browser()
  #Read grid, row and column from the object
  grids<-fData(x)[,"Grid"]
  rows<-fData(x)[,"Row"]
  cols<-fData(x)[,"Column"]
  nrReplicates<-dim(Data(x))[2]
  
  #Read the number of wells in each plate from the object
  nrWell <- prod(pdim(x))
  #Read the number of plates in the screen from the object
  nrPlates <- max(plate(x))
  #Compute the number of blocks from the objects position data
  nrGrids<-max(grids)
  #Read the number of channels from the object
  nrChannels<-dim(Data(x))[3] 
  
  #Extract the configured and annotated assay data into a data frame and add the plate and position data
  xAssayData<-Data(x)
  xAssay<-data.frame(Data(x))
  xAssay<-cbind(xAssay,plate=plate(x),grids=grids,blockRow=rows,blockCols=cols,order=order)
  
  #Create a boxplot of the raw pin block intensities of channel 1
  #boxplot(xAssay$X1.intensity~blocks, data=xAssay)
  
  #Normalize each channel and replicate in each block in each plate
    for(p in 1:nrPlates) {
      for(chr in 1:(nrChannels*nrReplicates)){
        for(b in 1:nrGrids){
          #browser()
              xAssay[,chr][xAssay$plate==p & xAssay$grids==b]<-
                pinBlock(xAssay[,chr][xAssay$plate==p & xAssay$grids==b],order, method,drawPlots)
            }
          }
        }
#      }
#    }
  
  #Put the normalized data back into the original array
  #This replaces each channel separately. Should be upgraded with a more compact method
  for (ch in 1:nrChannels){
    for(rep in 1:nrReplicates){
      xAssayData[,rep,ch]<-xAssay[,((ch-1)*nrReplicates+rep)]
    }
  }
  
#Put the array with the normalized data back into the new cellHTS object
  xpt<-x
  Data(xpt)<-xAssayData  
  xpt@state["normalized"]<-TRUE
  return(xpt)
}


#TODO: create a function that takes normalized and scored objects, creates a toptable
# and dedups it writing it to a file based on its name and returning a toptable dataframe

naiveScoreReplicates<-function(xsc){
  #Naive Replicate Scoring
  #Author: Mark Dane
  
  #Get the scores and gene ID's from the scored object
  xscData<-Data(xsc)
  xscfData<-fData(xsc)
  
  #Find all instances of the gene in the object and average their scores
  for(gene in xscfData$GeneID){
    same<-which(xscfData$GeneID %in% gene)
    score<-mean(xscData[same])
    xscData[same]<-score
  } 
  #Put the averaged scores back into the object
  Data(xsc)<-xscData
  
  #Return the scored object with average replicate scores
  return(xsc)
}

wellList<-function(xsc,thresh=c(-2.92,2.92)){
  #Make a dataframe from the assay and feature data
  xscd<-Data(xsc)
  xscf<-fData(xsc)
  out<-data.frame(xscf,xscd,seqPos=1:(length(xscd)/as.numeric(max(xscf$plate))))
  names(out)<-c('plate','well','finalWellAnno','siRNAID', 'GeneID','posName','score','position')
  
  #Eliminate genes with NA scores
  out<-out[!is.na(out$score),]
  #Identify which wells have duplicate genes and eliminate them
  dup<-duplicated(out$GeneID)
  outUnique<-out[which(!dup),]
  
  #Create a file for input to the scan^R
  #File format is: GeneID score well wellPos Plate finalWellAnno
  
  #Rearrange the output dataframe to create the file format order
  tempDF<-with(outUnique,data.frame(GeneID, score, well, wellPos=position,
                                    plate, finalWellAnno, stringsAsFactors=FALSE))
  
  #Sort the output by scores
  tempDF<-tempDF[order(tempDF$score, decreasing=TRUE),]
  
  #Apply the thresholds
  tempDF<-tempDF[(tempDF$score<=thresh[1]) | (tempDF$score>=thresh[2]),]
  
  #Write the dataframe to a tab delimited file with a header row
  write.table(tempDF,file='wellSelect.txt', sep="\t", row.names=FALSE, quote=FALSE)
  
  return(tempDF)
}

readBackwards<-function(x){
  #flip the columns in a dataframe up down and left right to make up for
  #an error in how the array-pro reads the data
  for(i in 1:dim(x)[2]){
    #reverse the order of the rows
    x[,i]<-x[dim(x)[1]:1,i]
  }
  return(x)
}

getANGrid<-function(grid, row, col, nrGridCols, nrRows, nrCols){
  #given vectors (grid, row, col) of well coordinates, the number of rows and columns in a grid (nrRows,nrcols),
  #and the number of grids in a CSMA row (nrGrids) this function returns the corresponding alhpanumeric notation
  plateRow<-as.integer((grid-1)/nrGridCols)*nrRows+row
  plateCol<-as.integer((grid-1)%%nrGridCols)*nrCols+col
  return(unname(unlist(getAlphaNumeric(plateRow,plateCol)[1])))
}

logObject<-function(x,base=2){
  #Returns the cellHTS2 object with the Data slot replaced with the log of the values
  #TODO: make this a robust log function to deal with non-positive values
  Data(x)<-log(Data(x), base=base)
  Data(x)[Data(x)<=-1e10]<-NA
  return(x)
}

selectColumn<-function(wd=getwd(),infile,outfile,sel=1){
  x<-read.table(file=paste(wd,infile,sep=""), header=TRUE, sep="\t",stringsAsFactors=FALSE)
  write.table(cbind(x[1],x[,sel],x[,dim(x)[2]]),
              file=paste(wd,outfile, sep=""),
              sep="\t", quote=FALSE,col.names = c("\t","intensities","position"))
}

orderDF<-function(df,blockGridCols=4){
  #Create an order array that converts from block, row, column space to subarray row and column space
  #This is used for GAL and aray pro files that are organized by block, row, columns
  #Return the dataframe in the correct space
  
  blockGridRows<-max(df$Grid)/blockGridCols
  blockRows<-max(df$Row)
  blockCols<-max(df$Column)
  blockSize<-blockRows*blockCols
  nrBlocks<-blockGridRows*blockGridCols
  nrWells<-nrBlocks*blockSize
  df<-cbind(df,Index=1:nrWells)
  
  #Create ordering pin block matrices.
  outDF<-df
  for(blockGridRow in 0:(blockGridRows-1)){
    blocks<-(((blockGridRow)*blockGridCols)+1):(((blockGridRow)*blockGridCols)+blockGridCols)
    selectBs<-df[df$Grid %in% blocks,]
    ord<-selectBs$Index[order(selectBs$Row,selectBs$Grid,selectBs$Column)]
    outDF[df$Grid %in% blocks,]<-df[ord,]
  }
  
return(outDF)
}

orderDFV2<-function(df,gridsPerRow=4){
  #Convert a dataframe from block, row, column space to subarray row and column space
  #This is used for GAL and aray pro files that are organized by block, row, columns
  #Return the dataframe in the correct space
  rowsPerGrid<-max(df$Row)
  colsPerGrid<-max(df$Column)
  df$arrayGridRow<-ceiling(df$Grid/gridsPerRow)
  df$arrayRow<-(df$arrayGridRow-1)*rowsPerGrid+df$Row
  df$arrayColumn<-((df$Grid-1)%%gridsPerRow)*colsPerGrid+df$Column
  df<-df[order(df$arrayRow,df$arrayColumn),]
  df<-subset(df,select=-c(arrayGridRow))
  df$PlateIndex<-1:(max(df$arrayColumn)*max(df$arrayRow))
return(df)
}

pseudoPrintOrder<-function(blockRows=10,blockCols=10,blockGridRows=12,blockGridCols=4){
  #return an array of a pseudo print order that is sequential within each pin grid
  blockSize<-blockRows*blockCols
  nrBlocks<-blockGridRows*blockGridCols
  nrWells<-nrBlocks*blockRows*blockCols
  
  #Create ordering pin block matrices
  blocks<-list()
  for(block in 1:nrBlocks){
    blocks[[block]]<-matrix(1:blockSize, ncol=blockCols, byrow=TRUE)
  }
  #Combine 4 pin block matrices to form a row
  orderMat<-cbind(blocks[[1]], blocks[[2]],
                  blocks[[3]], blocks[[4]])
  
  #Combine row groups into a matrix that matches the array
  for(blockGridRow in 2:blockGridRows){
    index<-(blockGridRow-1)*blockGridCols+1
    blockRow<-cbind(blocks[[index]], blocks[[index+1]],
                    blocks[[index+2]], blocks[[index+3]])
    orderMat<-rbind(orderMat,blockRow)
  }
  
  #Read the matrix out by row into an ordering array
  orderArray<-array(t(orderMat))
  return(orderArray)
}

orderArray<-function(positionRows,positionColumns,blockGridRows=12,blockGridCols=4){
  #Create an order array that converts from block, row, column space to subarray row and column space
  #This is used for GAL file annotations that are organized by block, row, columns
  #Using the output of this function to order a GAL file will put order it by rows and then columns
  #Caution on R's default of reading by column then row.
  #Usage restriction: this function assumes there are four pin grids across a subarray row
  #TODO: Make this general for rows that don't have 4 grids columns
  #TODO Generalize this for the number of plates
  
  blockRows<-max(positionRows)
  blockCols<-max(positionColumns)
  blockSize<-blockRows*blockCols
  nrBlocks<-blockGridRows*blockGridCols
  nrWells<-nrBlocks*blockSize
  
  #Create ordering pin block matrices. 
  blocks<-list()
  for(block in 1:nrBlocks){
    blocks[[block]]<-matrix(((block-1)*blockSize+1):((block-1)*blockSize+blockSize),
                            ncol=blockCols, byrow=TRUE)
  }
  #Combine 4 pin block matrices to form a row
  orderMat<-cbind(blocks[[1]], blocks[[2]],
                  blocks[[3]], blocks[[4]])
  
  #Combine row groups into a matrix that matches the array
  for(blockGridRow in 2:blockGridRows){
    index<-(blockGridRow-1)*blockGridCols+1
    blockRow<-cbind(blocks[[index]], blocks[[index+1]],
                    blocks[[index+2]], blocks[[index+3]])
    orderMat<-rbind(orderMat,blockRow)
  }
  
  #Read the matrix out by row into an ordering array
  orderArray<-array(t(orderMat))
  orderArray<-c(orderArray, nrWells+orderArray)
  return(orderArray)
}


logFileFUN<-function(x,thresh=2){
  #Return a dataframe for the log file that has the spots with  cell counting channel (3)
  #readings below a threshold. These will be flagged in the downstream analysis
  #Ouput format to include Plate, Well, Flag, Channel of the spots less than thresh
  
  #Get the data and feature data from the object. Returns as a 3D matrix
  #dim 1 is the rows
  #dim 2, ? and set to 1 for the debugging example
  #dim 3 has a dimension for each channel of data
  #TODO: Get this function to stop thowing warnings about the row names
  xData<-Data(x)
  xfData<-fData(x)
  #Create an empty data frame
  logDF<-data.frame()
  channel<-dim(Data(x))[3]
  #Try to use the channel data directly in xData to build the output dataframe
  #Cycle through all four replicates which are stored in the second dimension of xData
  for(rep in 1:dim(xData)[2]){
  #find spots in channel within the flagging range  
    #Don't attempt to make a dataframe if there are no values to be flagged
    if(length(which(xData[,rep,channel]<=thresh)) | sum(is.na(xData[,rep,channel]))){
      #Found at least one value in the flagging range so create a dataframe of all spots found
      foo<-data.frame(Plate=1,Sample=rep,
                      Well=xfData$well[which(xData[,rep,channel]<=thresh | is.na(xData[,rep,channel]))],
                      Flag="NA",Channel=channel,
                      stringsAsFactors=FALSE,row.names=NULL)
      #add the current spots to any already found
      logDF<-rbind(logDF,foo)
    }#End current replicate
  }#End for(rep in 1:dim(x)[2]){
  
  return(logDF)
}#end LogFileFUn

spatialVariation<-function(x){
  #Do a spatial variation analysis of the raw data of each channel in a cellHT object that has been 
  #normalized with the locfit method.
  #Get the locfit model of the raw data
  #######Hard coded this only works for cellHTS objects with one plate
  nrPlates<-length(x@plateData)
  nrChannels<-dim(Data(x))[3]
  nrReplicates<-dim(x@plateData[[1]])[2]
  
  #Do a spatial variation analysis of the raw data
  #Get the locfit model of the raw data
  pe<-plateEffects(x)
  
  ####Hard coded for 1 plate #####Get the residuals of the loc fit model from its median
  #This create a deviation value of the model from the plate median for every spot/well
  #This value should capture the magnitude of larger spatial effects while ignoring smaller ones
  #The nn parameter in locfit controls the sensitivity between large and small spatial effects
  spatialVarsReps<-sapply(X=1:nrReplicates,FUN=function(rep,x,pe,nrChannels){
    #Get the total variation of the model from it's median
    peChRes<-sapply(X=1:nrChannels,FUN=function(ch,pe,rep){
      #Get the variation of all channels in the current replicate
      #as a proportion of the median
      #Get a logical vector of the non-empty wells
      nonEmpty<-matrix(data=x@plateConf$Content!="empty",nrow=dim(x)[1],byrow=FALSE)
      #Remove the empty wells from pe and the calculation of Spatial Variation
       peRepCh<-pe$rowcol[nonEmpty[,rep],rep,ch]     
      return(sum(abs(peRepCh-median(peRepCh))/median(peRepCh)))
    },pe=pe,rep=rep)   
  },x=x,pe=pe,nrChannels=nrChannels)
  return(spatialVarsReps/dim(x)[1])
}

getRawExcelDF<-function(file, dataType="Net",...){
  #Read an excel file into a dataframe
  temp<-read.xls(file, verbose=FALSE, stringsAsFactors=FALSE,check.names=TRUE,...)
  #Find the dataType columns
  dataCols<-grep(pattern=dataType,x=names(temp))
  #Find the pin grid columns
  annotationCols<-which(names(temp) %in% c("Grid","Row","Column"))
  #Work around bug in read.xls that reads some numbers as characters
  df<-data.frame(unlist(apply(X=temp[,c(dataCols,annotationCols)],
                              MARGIN=2,FUN=as.numeric)))
  #Remove statistic rows after the data
  maxRow<-which(temp$X=="Maximum")
  df<-df[1:(df$Grid[maxRow]*df$Row[maxRow]*df$Column[maxRow]),]
  rm(temp)
  #Make meaningful short names for the columns of the channel number,
  #wavelength and dataType
  for(col in 1:length(dataCols)){
    #find the wavelength in the column name
    m <- regexpr(pattern="488|532|635", text=names(df[col]))
    #Get the wavelength from the column name
    wv<-regmatches(names(df[col]), m)
    names(df)[col]<-paste0("Ch",col,"_",wv,"_",dataType)
  }
  return(df)
}#End getRawexcelDF

# Gamma distribution approximation of p value for rank product.
#Author R.Eisinga et al./FEBS Letters 587 (2013) 677-682
righttailgamma = function(r,k,n) 1 - pgamma(-log(r/(n+1)^k),k,scale=1)

rbindDfs<-function(dfs){
  #rbind the replicates into a single dataframe and add the replicate column
  out<-dfs[[1]]
  for(df in dfs[-1]){
    out<-rbind(out,df)
  }
  out$replicate<-rep(x=1:length(dfs),each=dim(dfs[[1]])[1])
  return(out)
}#End rbindDFs

readReplicates<-function(dataFile,gal){
  #read the datafile and merge it with the gal annotations
  
  #Read the excel spreadsheet into a dataframe
  values<-read.xls(dataFile, verbose=FALSE, stringsAsFactors=FALSE)
  #Change the Grid column to numeric
  values$Grid<-as.numeric(values$Grid)
  #merge the gal and values dataframes on Grid, Row and Column
  values<-values[1:max(suppressWarnings(as.double(values$X)), na.rm=TRUE),]
  sheetDF<-merge(x=gal,y=values,by=c("Grid","Row","Column"))
  #Order the input dataframes by row then column
  sheetDF<-sheetDF[order(sheetDF$arrayRow,sheetDF$arrayColumn),]
  
  #Add a check on annotations from the gal file matching annotations from values
  
}#End readReplicates

readReplicatesAndGal<-function(dataFile){
  #read a datafile with  gal annotations
  
  #Read the excel spreadsheet into a dataframe
  values<-read.xls(dataFile, verbose=FALSE, stringsAsFactors=FALSE)
  #Change Block to Grid
  names(values)[which(names(values) %in% "Block")]="Grid"
  #Add plate and replicate columns
  values$plate=1
  #Add a check on annotations from the gal file matching annotations from values
  return(values)
}#End readReplicatesAndGal

simplifyNames<-function(df,pattern="Background|Net|Raw"){
  #Find columns with a key name and simplify it to it's type and wavelength
  colNames<-names(df)
  temp<-lapply(colNames,FUN=function(colName,pattern,df){
    m<-regexpr(pattern,colName)
    if(m[[1]]!=-1){
      #Found a column with one of the types of data in its name
      dt<-regmatches(colName,m)
      m<-regexpr(pattern="488|532|635|DAPI|594|647",colName)
      if(m[[1]]!=-1){
        #Found a column that has the right type of data and a wavelength
        wl<-regmatches(colName,m)
        #add a channel number prefix
        ######ch<-switch(EXPR=wl,"488"="Ch1","532"="Ch2","635"="Ch3","DAPI"="Ch1","594"="Ch2","647"="Ch3")
        #create the new name for the column/channel
        #####colName<-paste0(ch,dt,wl)
        colName<-paste0(dt,wl)
      }
    }
    return(colName)
  },pattern=pattern)
  names(df)<-temp
  return(df)
}#End simplifyNames

makeNumeric<-function(df,pattern="Ch|Grid|Net|Background|Raw|X"){
  #Find columns with a channel name
  columns<-grep(pattern,names(df))
  #Change each channel values to numeric
  for(c in columns){
    df[c]<-unlist(lapply(X=df[c],FUN=as.numeric))
  }
  return(df)
}#End makeNumeric

addControls<-function(df,pos=NULL,neg=NULL,empty=NULL){
  #add a configuration column to the dataframe
  #Wells are samples if they are not empty or a control
  df$controlStatus<-rep('sample', dim(df)[1])
  #Add any negative controls into this list
  selectNeg<-df$name %in% c('Negative', neg)
  selectNeg<-df$Name %in% c('Negative', neg)
  df$controlStatus[selectNeg]<-'neg'
  #Add any positive controls into this list
  selectPos<-df$name %in% c('Positive',pos)
  selectPos<-df$Name %in% c('Positive',pos)
  df$controlStatus[selectPos]<-'pos'
  selectEmpty<-df$name %in% c('Empty',empty)
  selectEmpty<-df$Name %in% c('Empty',empty)
  df$controlStatus[selectEmpty]<-'empty'
  return(df)
}#End addControls

filterOnCellCount<-function(df,countCh,thresh=0){
  #Replace the data values in all channels if the cell counts are below thresh or above the 99.9th percentile
  df<-Data(df)
  for (r in 1:dim(df)[2]){#go through all replicates
    for(row in 1:dim(df)[1]){#search in all rows of ch 1 of this replicate for low values
      #Force all low readings to NA values
      if(df[row,r,countCh]<=thresh |
           df[row,r,countCh] > quantile(df[,r,countCh],probs=.999, na.rm=TRUE)) df[row,r,]=NA
    }
  }  
  #Put the to cell count channel back into the cellHTS2 object
  Data(x)<-df
  return(x)
}#End filterLowCellCount

chtsAnnotate<-function(annDF,x){
  #annotatae the x object with the annDF dataframe
  tempFile<-tempfile()
  write.table(annDF,file=tempFile, quote=FALSE, sep="\t")
  x<-annotate(x, geneIDFile=tempFile)
  return(x)
}#End chtsAnnotate

chtsNormalizePinGrids<-function(x,pinGridNorm='none'){
  #Normalize all replicates and channels in x usingthe pinGridNorm method
  if(pinGridNorm=="none") {
    xpt<-x
  } else {
    xpt<-pinTipNormalize(x,method=pinGridNorm)
  }
}#End chtsNormalizePinGrids

chtsNormalizePlates<-function(x,normMethod='none',logged=FALSE,refChannel=NULL,...){
  #Wrapper to add reference channel method to normalize plates
  if (normMethod!='none'){
    if (normMethod=='refChannel') x<-refChannelNormalize(x, refChannel=refChannel,...)
    else {
      if(logged){
        x<- suppressWarnings(normalizePlates(x,scale="additive", method=normMethod, ...))
      }else {
        x <- suppressWarnings(normalizePlates(x,scale="multiplicative", method=normMethod,...))
      }
    }
  }
  return(x)
}#End chtsNormalizePlates

trimSummaryStats<-function(df){
  #Remove the summary statistics from the end of ProArray dataframe
  if(length(which(df$X=="Mean value"))!=0){
    tmp<-df[1:(which(df$X=="Mean value")-1),]
  } else {
    tmp<-df
  }
  return(tmp)
}#End trimSummaryStats

getAnnData<-function(x){
  #extract the data, plate, grid(Grid), row and column from a cellHTS2 object and return it in a dataframe
  #There is one column for each channel. Replicate values denoted by the replicate column
  nrChannels<-dim(Data(x))['Channels']
  nrReplicates<-dim(Data(x))['Samples']
  channels<-data.frame(Data(x))
  annotations<-fData(x)
  #Find the replicate 1 channels by their names starting with 'X1.'
  pattern<-'X1.'
  repChannels <-grep(pattern,names(channels))
  #Create a dataframe with the replicate one channels, the annotations and replicate =1
  xout<-cbind(channels[repChannels],annotations,replicate=1)
  #Strip off the X1. from the channel names
  names(xout)[1:nrChannels]<-gsub("X[0-9].", "",names(channels))[repChannels]
  #Repeat the same process for the rest of the replicates, and add the to the bottom of the dataframe
  if (nrReplicates>1){
    for(r in 2:nrReplicates){
      pattern<-paste('X',r,sep="")
      repChannels <-grep(pattern,names(channels))
      xrep<-cbind(channels[repChannels],annotations,replicate=r)
      names(xrep)[1:nrChannels]<-gsub("X[0-9].", "",names(channels))[repChannels]
      xout<-rbind(xout,xrep)
    }
  }
  xout$controlStatus<-factor(xout$controlStatus, levels=c('sample','pos','neg','other'), order=TRUE)
  return(xout)
}

chtsPinBoxPlots<-function(x, main="",col='grey',...){
  #Create boxplots of the pin tip values from a cellHTS2 object
  xdf<-getAnnData(x)
  #Read the number of wells in each plate from the object
  nrWell <- prod(pdim(x))
  #Read the number of plates in the screen from the object
  nrPlates <- max(plate(x))
  #Compute the number of Grids from the objects position data
  nrGrids<-max(xdf$Grid)
  #Read the number of channels from the object
  nrChannels<-dim(Data(x))[3] 
  #Read the number of replicates
  nrReplicates<-dim(Data(x))['Samples']
  for(p in 1:nrPlates) {
    for(r in 1:nrReplicates){
      for(ch in 1:nrChannels){
        boxplot(xdf[,ch][xdf$replicate==r]~xdf$Grid[xdf$replicate==r],data=xdf, groups=xdf$controlStatus,
                col=col, aspect = "xy",
                main=main,
                #paste(main,'\n','Replicate', r, 'Channel',ch, 'Layout',p),
                xlab='Pin Grid',ylab='Channel Value',...)
      }
    }
  }
}


chtsPinQQs<-function(x,main=""){
  #Create qqplots of the pin grids of a cellHTS2 object with different colors for the sample types. 
  library(lattice)
  xdf<-getAnnData(x)
  #Read the number of wells in each plate from the object
  nrWell <- prod(pdim(x))
  #Read the number of plates in the screen from the object
  nrPlates <- max(plate(x))
  #Compute the number of Grids from the objects position data
  nrGrids<-max(xdf$Grid)
  #Read the number of channels from the object
  nrChannels<-dim(Data(x))[3]
  for(p in 1:nrPlates) {
    #Plates are not implemented yet(?)
    for(r in 1:nrReplicates){
      for(ch in 1:nrChannels){
        plotqq<-qqmath(~ xdf[,ch][xdf$replicate==r] | Grid[xdf$replicate==r] , data=xdf, groups=controlStatus,
                       aspect = "1",pch=c(21,21,21),col=c('grey','black','black','black'),
                       fill=c('black','red','blue','magenta2'),
                       main=paste(main,'Pin Tip','\n',
                                  'Replicate',r,'Channel',ch,'Plate',p), strip=strip.custom(strip.names = c(FALSE, FALSE)),
                       ylab='intensity',  
                       prepanel = function(x, ...) {
                         list(xlim = range(qnorm(ppoints(length(x)))))
                       },layout=c(4,12),
                       panel = function(x, ...) {
                         qx <- qnorm(ppoints(length(x)))[rank(x)]
                         panel.xyplot(x = qx, y = x,..., aspect="1")
                         panel.abline(coef(lm(x~qx)), col='red', aspect="1")
                       })
        print(plotqq)
      }
    }
  }
}

chtsPlateQQs<-function(x,main=""){
  #Create qqplots of the pin grids of a cellHTS2 object with different colors for the sample types. 
  library(lattice)
  xdf<-getAnnData(x)
  #Read the number of wells in each plate from the object
  nrWell <- prod(pdim(x))
  #Read the number of plates in the screen from the object
  nrPlates <- max(plate(x))
  #Compute the number of Grids from the objects position data
  nrGrids<-max(xdf$Grid)
  #Read the number of channels from the object
  nrChannels<-dim(Data(x))[3]
  for(p in 1:nrPlates) {
    #Plates are not implemented yet(?)
    for(r in 1:nrReplicates){
      for(ch in 1:nrChannels){
        plotqq<-qqmath(~ xdf[,ch][xdf$replicate==r] , data=xdf, groups=controlStatus,
                       aspect = "xy",pch=c(21,21,21),col=c('black','black','black','black'),
                       fill=c('black','red','blue','magenta2'),
                       main=paste(main,'Plate Values', '\n',
                                  'Replicate ',r,'Channel',ch,'Plate',p), strip=strip.custom(strip.names = c(FALSE, FALSE)),
                       ylab='intensity',  
                       prepanel = function(x, ...) {
                         list(xlim = range(qnorm(ppoints(length(x)))))
                       },layout=c(1,1),
                       panel = function(x, ...) {
                         qx <- qnorm(ppoints(length(x)))[rank(x)]
                         panel.xyplot(x = qx, y = x,...)
                         panel.abline(coef(lm(x~qx)), col='red')
                       })
        print(plotqq)
      }
    }
  }
}

chtsPinHists<-function(x,main=""){
  #Create histograms of the pin grids of a cellHTS2 object
  #TODO add a rug with different colors for the sample types
  library(lattice)
  xdf<-getAnnData(x)
  #Read the number of wells in each plate from the object
  nrWell <- prod(pdim(x))
  #Read the number of plates in the screen from the object
  nrPlates <- max(plate(x))
  #Compute the number of Grids from the objects position data
  nrGrids<-max(xdf$Grid)
  #Read the number of channels from the object
  nrChannels<-dim(Data(x))[3]
  for(p in 1:nrPlates) {
    #Plates are not implemented yet(?)
    for(r in 1:nrReplicates){
      for(ch in 1:nrChannels){
        rawhist<-histogram(~xdf[,ch][xdf$replicate==r] | Grid[xdf$replicate==r], data=xdf,
                           groups=controlStatus,
                           layout = c(4,2), 
                           xlab = "Intensity",
                           main=paste(main,'Pin Tip Values','\n',
                                      'Replicate', r, 'Channel',ch,'Plate',p),
                           panel = function(x, ...) {
                             panel.histogram(x, col='grey', ...)
                             panel.rug(x,col='black', ...)
                           }  )
        print(rawhist)
      }
    }
  }
}

chtsPlateHists<-function(x,main=""){
  #Create histograms of the pin grids of a cellHTS2 object
  #TODO add a rug with different colors for the sample types
  library(lattice)
  xdf<-getAnnData(x)
  #Read the number of wells in each plate from the object
  nrWell <- prod(pdim(x))
  #Read the number of plates in the screen from the object
  nrPlates <- max(plate(x))
  #Compute the number of Grids from the objects position data
  nrGrids<-max(xdf$Grid)
  #Read the number of channels from the object
  nrChannels<-dim(Data(x))[3]
  for(p in 1:nrPlates) {
    #Plates are not implemented yet(?)
    for(r in 1:nrReplicates){
      for(ch in 1:nrChannels){
        rawhist<-histogram(~xdf[,ch][xdf$replicate==r] , data=xdf,
                           groups=controlStatus,
                           layout = c(1,1), 
                           xlab = "Intensity",
                           main=paste(main,'Plate Values','\n',
                                      'Replicate', r, 'Channel',ch,'Plate',p),
                           panel = function(x, ...) {
                             panel.histogram(x, col='grey', ...)
                             panel.rug(x,col='black', ...)
                           }  )
        print(rawhist)
      }
    }
  }
}

chtsPinXY<-function(x, main=""){
  #Create scatterplots of the pin grids of a cellHTS2 object
  library(lattice)
  #Create scatterplots of the raw pin grid values
  xdf<-getAnnData(x)
  #Read the number of wells in each plate from the object
  nrWell <- prod(pdim(x))
  #Read the number of plates in the screen from the object
  nrPlates <- max(plate(x))
  #Compute the number of Grids from the objects position data
  nrGrids<-max(xdf$Grid)
  #Read the number of channels from the object
  nrChannels<-dim(Data(x))[3]
  GridRows=max(xdf$Row)
  GridCols<-max(xdf$Column)
  #Assume there are 4 Grid columns in a subarray
  GridGridCols=4
  #First create a print order array
  temp<-pseudoPrintOrder(blockCols=GridCols, blockRows=GridRows,
                         blockGridRows=nrGrids/GridGridCols, blockGridCols=GridGridCols)
  for(p in 1:nrPlates) {
    for(r in 1:nrReplicates){
      for(ch in 1:nrChannels){
        XY<-xyplot(xdf[,ch][xdf$replicate==r] ~ temp | Grid[xdf$replicate==r], data=xdf,
                   groups=controlStatus,pch=c(21,21,21,21),col=c('grey','red','blue','magenta2'),
                   fill=c('black','red','blue','magenta2'),layout = c(2,2), 
                   xlab = "Order", 
                   ylab='Intensity',
                   main=paste(main,'Pin Tip Values','\n',
                              'Replicate', r, 'Channel',ch,'Plate',p),
                   panel = function(x, y,...) {
                     panel.xyplot(x, y,...)
                     panel.abline(median(y),0, col='blue')
                     panel.loess(x,y, col='red')
                     
                   })
        print(XY)
      }
    }
  }
}

chtsPlateXY<-function(x, main=""){
  #Create scatterplots of the pin grids of a cellHTS2 object
  library(lattice)
  #Create scatterplots of the raw pin grid values
  xdf<-getAnnData(x)
  #Read the number of wells in each plate from the object
  nrWell <- prod(pdim(x))
  #Read the number of plates in the screen from the object
  nrPlates <- max(plate(x))
  #Compute the number of Grids from the objects position data
  nrGrids<-max(xdf$Grid)
  #Read the number of channels from the object
  nrChannels<-dim(Data(x))[3]
  GridRows=max(xdf$Row)
  GridCols<-max(xdf$Column)
  #Assume there are 4 Grid columns in a subarray
  GridGridCols=4
  #First create a print order array
  temp<-pseudoPrintOrder(blockCols=GridCols, blockRows=GridRows,
                         blockGridRows=nrGrids/GridGridCols, blockGridCols=GridGridCols)
  for(p in 1:nrPlates) {
    for(r in 1:nrReplicates){
      for(ch in 1:nrChannels){
        XY<-xyplot(xdf[,ch][xdf$replicate==r] ~ temp , data=xdf,
                   groups=controlStatus,pch=c(21,21,21,21),col=c('grey','red','blue','magenta2'),
                   fill=c('black','red','blue','magenta2'),layout = c(1,1), 
                   xlab = "Order", 
                   ylab='Intensity',
                   main=paste(main,'Plate Values','\n',
                              'Replicate', r, 'Channel',ch,'Plate',p),
                   panel = function(x, y,...) {
                     panel.xyplot(x, y,...)
                     panel.abline(a=median(y),b=0, col='blue')
                     panel.loess(x,y, col='red')
                     
                   })
        print(XY)
      }
    }
  }
}


chtsPlateByPinXY<-function(x, main=""){
  #Create scatterplots of the pin grids of a cellHTS2 object
  library(lattice)
  #Create scatterplots of the raw pin grid values
  xdf<-getAnnData(x)
  for(p in 1:nrPlates) {
    for(r in 1:nrReplicates){
      for(ch in 1:nrChannels){
        XY<-xyplot(xdf[,ch][xdf$replicate==r] ~ Grid , data=xdf,
                   groups=controlStatus,pch=c(21,21,21,21),col=c('grey','red','blue','magenta2'),
                   fill=c('black','red','blue','magenta2'),layout = c(1,1), 
                   xlab = "Order", 
                   ylab='Intensity',
                   main=paste(main,'Plate Values','\n',
                              'Replicate', r, 'Channel',ch,'Plate',p),
                   panel = function(x, y,...) {
                     panel.xyplot(x, y,...)
                   })
        
        print(XY)
      }
    }
  }
}

########## Rank Product of Replicates #########

chtsRankReplicates<-function(x,channels){
  #Rank each replicate on a per channel basis, calculate the rank products
  #across the replicates then store the p values of the rank products in the
  #data slot, overwriting the normalized data
  
  #Get the assay data
  xd<-Data(x)
  
  #Check the channel names and get the indices
  chs<-which(channelNames(x) %in% channels)
  if(length(chs)==0)stop("no channels matched the requested ones.")
  
  #Get the number of replicates and check its more than 1
  reps<-as.numeric(unlist(dimnames(xd)[2]))
  if(length(reps)<=1)stop("Rank products require two or more replicates")
  
  #Check the state of the object
  if(!all(state(x)[c(1,2,4)]))
    stop("The cellHTS2 object must be configured, normalized and annotated")
  ranks<-matrix(nrow=dim(xd)[1],ncol=length(reps),
                dimnames=list(NULL,paste0("ranks_r",colnames(xd))))
  #order both replicates by their channel values
  for(r in reps){
    #The number one rank will have the highest value
    ranks[,r]<-(1+dim(xd)[1])-rank(xd[,r,chs])
  }
  
  #Calculate the rank products
  rp<-apply(X=ranks,MARGIN=1,FUN=prod)
  
  #Get the rank product p value estimates using the gamma distribution
  #This can be replaced with the exact p-values if needed
  #rpp<-sapply(X=rp,FUN=righttailgamma,k=length(reps),n=dim(xd)[1])
  k=length(reps)
  n=dim(xd)[1]
  #Debug Return the lower of the right tail or left tail p value
  rpp<-pgamma(-log(rp/(n+1)^k),k,scale=1)
  
  #Add the rpp scores into the object by replacing the normalized data
  x<-summarizeReplicates(x)
  Data(x)<-array(data=rpp,dim=c(length(rpp),1,1))
  
  return(x)
}
########## End Rank Product of Replicates #######
