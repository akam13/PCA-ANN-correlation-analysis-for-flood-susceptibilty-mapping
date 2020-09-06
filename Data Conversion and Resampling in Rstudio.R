##References
##http://neondataskills.org/R/Raster-Data-In-R/
##http://r-sig-geo.2731867.n2.nabble.com/How-I-make-2-rasters-with-equal-extents-td7584918.html
##https://geoscripting-wur.github.io/IntroToRaster/


     #################################      
##################  Omar F. Althuwaynee  ###################
################  PhD, GIS and Geomatics Engineering########
     #######          Mungu Betrand A    #############
     ###### Phd Candidate,Surveying and Mapping ######
     #################################################

     
 ###### PART 1 : Project General Settings #####
                  ##################################

#######################################################
###1 collecting all the results in one single folder###
#######################################################
  #create new two folders, 1)Working directory 2)Personallibrary
dir.create("C:/Users/user/Desktop/WorkingDIR") # change the / only and add ""
dir.create("C:/Users/user/Desktop/WorkingDIR/MyLibrary")  
#Set working directory (work space)
setwd("C:/Users/user/Desktop/WorkingDIR")
getwd()
 ##set the personal library
.libPaths("C:/Users/user/Desktop/WorkingDIR/MyLibrary")
.libPaths() ## Press tab for options










###################################
###2 Install required packages (Raster, rgeos, rgdal, maptools###
##################################
#To work with rasters in R, we need two key libraries, 
#sp and raster.rgdal will allow us to export rasters to geotiff format.
# install single packages
install.packages("raster")
#read package
library(raster)


# install and read maultiple packages code
packages <- c("rgeos", "rgdal", "maptools")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(packages)
##########################################################












###################################
###3 load raster in an R object###
##################################

#load raster in an R object called 'elevation'
Elevation = raster("MyData/Elevation.tif")  
list.files("MyData")
# look at the raster attributes. 
plot(Elevation, main="Elevation map")
#view coordinate reference system
Elevation@crs
#we can look at the distribution of values in the raster too
hist(Elevation, main="Distribution of elevation values", 
     col= "yellow", maxpixels=2000000)
#multiple each pixel in the raster by 5
Elevation_by5 <- Elevation * 5


#load all rasters from working directory 
Elevation = raster("MyData/Elevation.tif")  
Ndvi= raster("MyData/Ndvi.tif")
dfr = raster("MyData/dfr.tif")
slope= raster("MyData/slope.tif") 
twi= raster("MyData/twi.tif") 
rain = raster("MyData/rain.tif")
soil = raster("MyData/soil.tif")
spi = raster("MyData/spi.tif")
lulc = raster("MyData/lulc.tif")
class = raster("MyData/class.tif")

plot(Rainn)
list.files("MyData") 

        ########################### End of Par1  ########################



                       ###########   DATA PREPRATION   #########
                               #############################
#Prepare independents factors (Slope, elevation, curvature, NDVI) 
# dependent factor (Landslides locations>> Raster structure format, in sake of data homogeneity)
#load raster in an R object called 'elevation'
list.files("MyData") # see the files names
Slope <- raster("MyData/Slope.tif")  
# look at the raster attributes. 
plot(Slope, main="Slope map")
# Reclassify to build categorial map, hard classification
m= c(0,7,1,10,30,2)
mat=matrix(m,ncol=2, byrow = TRUE)
slopeCat=reclassify(Slope, mat)
plot(slopeCat)
slopeCat
writeRaster(slopeCat,"slopeCat.tif", overwrite=TRUE)
Slope



# load all the data
Elevation = raster("MyData/Elevation.tif")  
Ndvi= raster("MyData/Ndvi.tif")
dfr = raster("MyData/dfr.tif")
slope= raster("MyData/slope.tif") 
twi= raster("MyData/twi.tif") 
rain = raster("MyData/rain.tif")
soil = raster("MyData/soil.tif")
spi = raster("MyData/spi.tif")
lulc = raster("MyData/lulc.tif")
class= raster("MyData/class.tif") 
list.files("MyData")

# if you have diffrent extent, then try to Resample them 
class_re <- resample(class,Elevation, resample='bilinear') 
Ndvi_re<- resample(Ndvi,Elevation, resample='bilinear')
dfr_re<- resample(dfr,Elevation, resample='bilinear')
slope_re<- resample(slope,Elevation, resample='bilinear')
twi_re<- resample(twi,Elevation, resample='bilinear')
rain_re<- resample(rain,Elevation, resample='bilinear')
soil_re<- resample(soil,Elevation, resample='bilinear')
spi_re<- resample(spi,Elevation, resample='bilinear')
lulc_re<- resample(lulc,Elevation, resample='bilinear')


# create resampled Rasters
dir.create("C:/Users/user/Desktop/WorkingDIR/Resampled data")  
writeRaster(class_re,"Resampled data/class.tif", overwrite=TRUE)
writeRaster(slope_re,"Resampled data/slope.tif", overwrite=TRUE)
writeRaster(soil_re,"Resampled data/soil.tif", overwrite=TRUE) # working!!
writeRaster(Ndvi_re,"Resampled data/Ndvi.tif", overwrite=TRUE)
writeRaster(dfr_re,"Resampled data/dfr.tif", overwrite=TRUE)
writeRaster(spi_re,"Resampled data/spi.tif", overwrite=TRUE)
writeRaster(twi_re,"Resampled data/twi.tif", overwrite=TRUE)
writeRaster(rain_re,"Resampled data/rain.tif", overwrite=TRUE)
writeRaster(lulc_re,"Resampled data/lulc.tif", overwrite=TRUE)
writeRaster(Elevation,"Resampled data/Elevation.tif", overwrite=TRUE)












#############   preparing data for Statisitcal analysis #########
              #############################

## Let us call the resampled rasters files
list.files("C:/Users/user/Desktop/WorkingDIR/Resampled data")
class=raster("Resampled data/class.tif")
Elevation=raster("Resampled data/Elevation.tif")
Ndvi=raster("Resampled data/Ndvi.tif") 
dfr=raster("Resampled data/dfr.tif")
slope=raster("Resampled data/slope.tif")
rain=raster("Resampled data/rain.tif")
soil=raster("Resampled data/soil.tif")
spi=raster("Resampled data/spi.tif")
lulc=raster("Resampled data/lulc.tif")
twi=raster("Resampled data/twi.tif")


hist(class)  #Plot histogram to see the shape of data distribution 

plot(Slope)
plot(class, add=TRUE)
## stack multiple raster files

Stack_List= list.files(path = "C:/Users/user/Desktop/WD/Resampled data",pattern = "tif$", full.names = TRUE)
Rasters=stack(Stack_List)
names(Rasters)

head(Rasters)   ## our stack variable, we see NA values, #### 
                ##that needto be removed or replaced by -9999###

value_table=getValues(Rasters)
head(value_table, n=6)

value_table=na.omit(value_table)
value_table=as.data.frame(value_table)
head(value_table, n=6)

#Export your values to text file, thats for double check or 
#use it in microsoft Excel for further calculations##
write.table(value_table, "mydata.txt", sep="\t") # export to text file
str(value_table)

##### DONE..! ####
## the data exported to Excel can be used for any kind of statistical or machine learning anlysis in our cas we use the exported data for flood susceptibilty mapping with ANN and PCA##







