## Luke Ozsanlav-Harris
## 21/04/2022
## Trim landcover images to just the extent of Islay as easier to work with smaller files

## Packages required
library(tidyverse)
library(data.table)
library(raster)



##
#### 1. Read in all the raster layers ####
##


## read in the raster with the habitat data
LandCov2015 <- raster("Landcover Data/RAW landcover data/lcm2015gb25m.tif")
LandCov2017 <- raster("Landcover Data/RAW landcover data/gb2017lcm25m.tif")
LandCov2018 <- raster("Landcover Data/RAW landcover data/gb2018lcm25m.tif")
LandCov2019 <- raster("Landcover Data/RAW landcover data/gb2019lcm25m.tif")
LandCov2020 <- raster("Landcover Data/RAW landcover data/gb2020lcm25m.tif")

## plot one of the images
plot(LandCov2015)



##
#### 2. Crop the Rasters to just the extent of Islay and west end of Jura ####
##

## Set the extent required
Islay_ext <- extent(114500, 160000, 639000,680000) # This extent applies for the CRS of the CEH land use data

## crop each of the rasters
LandCov2015_Islay <- crop(LandCov2015, Islay_ext)
LandCov2017_Islay <- crop(LandCov2017, Islay_ext)
LandCov2018_Islay <- crop(LandCov2018, Islay_ext)
LandCov2019_Islay <- crop(LandCov2019, Islay_ext)
LandCov2020_Islay <- crop(LandCov2020, Islay_ext)

## plot a couple of the rasters
par(mfrow = c(1, 2))
plot(LandCov2015_Islay); plot(LandCov2020_Islay)


## create a raster to look at which areas have changed habitat
## create an empty raster first
Change_raster <- LandCov2015_Islay
values(Change_raster) <- NA
values(Change_raster) <- ifelse(values(LandCov2019_Islay)==values(LandCov2020_Islay), 0, 1)

##plot the results
par(mfrow = c(1, 1))
plot(Change_raster)



##
#### 3. Save the cropped Rasters ####
##


writeRaster(LandCov2015_Islay, file = "Landcover Data/Islay landcover data/LandCov2015_Islay", overwrite=TRUE)
writeRaster(LandCov2017_Islay, file = "Landcover Data/Islay landcover data/LandCov2017_Islay", overwrite=TRUE)
writeRaster(LandCov2018_Islay, file = "Landcover Data/Islay landcover data/LandCov2018_Islay", overwrite=TRUE)
writeRaster(LandCov2019_Islay, file = "Landcover Data/Islay landcover data/LandCov2019_Islay", overwrite=TRUE)
writeRaster(LandCov2020_Islay, file = "Landcover Data/Islay landcover data/LandCov2020_Islay", overwrite=TRUE)





