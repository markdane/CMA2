#This is a parameter packet that creates and analyzes simulated data that as no variances


#Name the experiment
experimentName <- "Simulated Data with pin and plate variances"

#Create seed(s) for the random number generators
rGenSeed<-1234

#Size the pin head which creates the grids
#the pin structure determines tha plate layout. Each pin creates one grid.
plateGridCols<-4
plateGridRows<-12
#Size each grid
gridRows<-10
gridCols<-9
nrPins<-plateGridCols*plateGridRows

#Calculate the plate sizes for later use
plateRows<-plateGridRows*gridRows
plateCols<-plateGridCols*gridCols
nrWells<-plateRows*plateCols

pinToPinSD<-1e-4
scale<-sqrt(pinToPinSD)
shape<-1/(scale)
set.seed(rGenSeed)
pinMeans<-rgamma(nrPins, shape=shape,scale=scale)
pinSD<-1e-4
pinSDs<-pinSD*pinMeans


#Configure the plates, channels and replicates
nrPlates <- 1
plateToPlateSD<-1e-10
scale<-sqrt(plateToPlateSD)
shape<-1/scale
set.seed(rGenSeed)
plateMeans<-rgamma(nrPlates, shape=shape, scale=scale)
plateSD<-1e-10
plateSDs<-plateSD*plateMeans

#Set the number of replicate plates
nrReplicates <- 2
repToRepSD<-1e-4
scale<-sqrt(repToRepSD)
shape<-1/scale
set.seed(rGenSeed)
replicateMeans<-rgamma(nrReplicates,shape=shape, scale=scale)
replicateSD<-1e-1
replicateSDs<-replicateSD*replicateMeans

#Determine whether to log the raw values before all processing
logInts<-TRUE
if(logInts) {plateNormScale<-'additive'
}else {plateNormScale<-'multiplicative'
}
#plate and pin normalization methods
normMethod<-'locfit'
nn<-0.004
nnList<-10^seq(from=-1.6,to=-2.5, by=-.1)
pinTipMethod<-'none'

#Use the following value to get identical replicates
#replicateMeans<-rep(1,nrReplicates)
#replicateSDs<-rep(0,nrReplicates)
nrChannels<-1
channelMeans<-c(5000,7000, 10000)
channelSD<-0.1
channelSDs<-channelSD*channelMeans

#Set the genes of interest percentage
pctGenesOfInterest<-0.01

#Select the description file
descripFile<-'Simulated Description.txt'

#Create the factor to knockdown the genes of interest
# GOIFactor<-0.6
# GOISD<-0.05
#Test roc by having no difference in the goi knockdown
GOIFactor<-0.3
GOISD<-0.05

#Create the factor to knockdown the positive controls
posFactor<-0.6
posSD<-0.05

#Create the ellipsoid perturbations
#xc is the x (column) location of the center
#yc is the y (row) location of the center
#ax is the radius along the x (column) axis
#ay is the radius along the y (row) axis
#p is the highest amount of perturbation
#put an ellipsoid of positive perturbations at the center of the array
e1<-c(plateCols/2,plateRows/2,plateCols/4,plateRows/4,1)
#put an ellipsoid of negative perturbations at the lower, left corner of the array
e2<-c(0,plateRows,plateCols/4,plateRows/4,-1)
#Make a list of the ellipsoid definitions
ellipList<-list(e1,e2)

#Add line pertrbations
l1<-c(plateRows*.6,1,plateRows*.8,plateCols*.2,1,1)
l2<-c(plateRows*.8,plateCols*.2,plateRows*.9,plateCols*.2,1,.5)
l3<-c(plateRows*.9,plateCols*.2,plateRows,plateCols,1,1)
lineList<-(list(l1,l2,l3))

