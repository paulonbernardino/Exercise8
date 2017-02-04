###### WUR Geo-scripting course
### February 3rd 2017
### Paulo Bernardino
### Exercise 8

# Load libraries
library(sp)
library(raster)
library(randomForest)
library(rgdal)

# Load data
load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")

## clean water, clouds and cloud shadow from VCF data 
vcfGewata[vcfGewata>100]<-NA
hist(vcfGewata)

## Build a brick containing all data
gwt <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7, vcfGewata)
names(gwt) <- c("band1", "band2", "band3", "band4", "band5", "band7", "vcf")

## Extract all data to a data.frame
gwt.df <- as.data.frame(getValues(gwt)) 
head(gwt.df)
str(gwt.df)

## Assign bands and VCF values to vectors 
band1<-gwt.df$band1
band2<-gwt.df$band2
band3<-gwt.df$band3
band4<-gwt.df$band4
band5<-gwt.df$band5
band7<-gwt.df$band7
vcf<-gwt.df$vcf

## Visualizing the data and the relationship between bands and VCF tree cover
plot(gwt) 
summary(gwt)
hist(gwt)
pairs(gwt)

## Creat a lm model
lm1<-lm(vcf~band1+band2+band3+band4+band5+band7)
summary(lm1)
## From the results, we can see that band 7 is of little importance 
## in predicting VCF. For instance, if we remove band 7, the resulting
## R-squared would still be the same:
summary(lm(vcf~band1+band2+band3+band4+band5)) 

## Predict tree cover from model
predTC<- predict(gwt,model=lm1,na.rm=TRUE)
par(mfrow=c(1,1))
predTC[predTC<0]<-NA
plot(predTC,xlim=c(808755,855195),ylim=c(817635,852945))
str(predTC)

###### randomForest model ###

str(gwt)

## Load the calibration points
samples<-readOGR("data/calib.gwt.kml",layer="calib.gwt")
class(samples)
## Re-project spatialPointsDataFrame
samples2<-spTransform(samples,CRS(proj4string(gwt)))
samples2@coords<-coordinates(samples2)[,-3]
## Extract the surface reflectance
calib<- extract(gwt,samples2,df=TRUE)
## Combine the created DF to the description columm of the calibration dataset
calib2<- cbind(samples2$Description,calib)
calib2<-na.omit(calib2)
str(calib2)

## Calibrate a random forest model using the extracted data frame
rFm<- randomForest(vcf~band1+band2+band3+band4+band5+band7,data=calib2)

## Use the model to predict tree cover
tcMap<-predict(gwt,model=rFm)
plot(tcMap)

## Function to compute the RMSE
RMSE<-function(x,y) {
  sqrt(mean((y-x$vcf)^2,na.rm=TRUE))
}

## Compute the RMSE for the lm model
predTC.df<-as.data.frame(getValues(predTC))
RMSE(gwt.df,predTC.df)

## Compute the RMSE for the random forest model
tcMap.df<-as.data.frame(getValues(tcMap))
RMSE(gwt.df,tcMap.df)

## Compute RMSE for each class
## Make an NA-value raster based on the LC raster attributes
formask <- setValues(raster(samples2), NA)
## Assign 1 to formask to all cells corresponding to the forest class
formask[samples2$Description=="forest"] <- 1
## Create a grass mask
grassmask <- setValues(raster(samples2), NA)
grassmask[samples2$Description=="grass"] <- 1
##Create a crop mask
cropmask <- setValues(raster(samples2), NA)
cropmask[samples2$Description=="crop"] <- 1

## Get values for each class
tcfor<-crop(tcMap,formask)
tcfor.df<-as.data.frame(getValues(tcfor))
tcgrass<-crop(tcMap,grassmask)
tcgrass.df<-as.data.frame(getValues(tcgrass))
tccrop<-crop(tcMap,cropmask)
tccrop.df<-as.data.frame(getValues(tccrop))

## Compute RMSE for each class
RMSE(gwt.df,tcfor.df)
RMSE(gwt.df,tcgrass.df)
RMSE(gwt.df,tccrop.df)