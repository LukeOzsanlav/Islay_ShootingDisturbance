## Luke Ozsanlav-Harris
## 22/04/2022
## Combine GWfG and GBG tracking data and label with landcover classes


## Packages required
library(tidyverse)
library(lubridate)
library(data.table)
library(zoo)
library(geosphere)
library(sf)
library(sfheaders)
library(raster)



##
#### 1. Read in Datasets ####
##

##                                           ##
#### 1.1  Read in and clean Ecotone GPS data ####
##                                           ##


## set working directory and read in ecotone data
eco_data <- readRDS("tracking data/Eco_data_cleaned.RDS")

#### Data Cleaning 
## extracting just the GPS data
eco_gps <- filter(eco_data, is.na(Acc_Z) ==T)
rm(eco_data)

## parse timestamp
eco_gps$UTC_datetime <- ymd_hms(eco_gps$UTC_datetime)

## Delete observations where missing lat or long or a timestamp. 
ind <- complete.cases(eco_gps[, c("Latitude", "Longitude", "UTC_datetime")])
table(ind) # number of locations with missing coordinates or timestamp (if any)
eco_gps <- eco_gps %>% filter(ind) # removing rows with missing data (if any)

##Check for duplicated observations (ones with same lat, long, timestamp, and device ID)
ind2 <- eco_gps %>% 
        dplyr::select(UTC_datetime, Longitude, Latitude, Tag_ID) %>%
        duplicated()
sum(ind2) 
eco_gps <- eco_gps %>% filter(!ind2)

## Streamline the dataset for a first pass, first just have subset of the columns
gps <- subset(eco_gps, select = c("Tag_ID" , "UTC_datetime", "Latitude", "Longitude"))

## Now filter for the winter
## create a year day column first
gps$year_day <- yday(gps$UTC_datetime)
winter_eco_gps <- filter(gps, year_day >= 305 | year_day <= 92)
table(winter_eco_gps$Tag_ID)







##                                            ##
#### 1.2  Read in and clean Ornitela GPS data ####
##                                            ##

## Data Cleaning, extracting just the GPS data
orn_gps <- readRDS("tracking data/Ornitela_GPSClean.RDS")

## read in meta data so i can just keep the Scotland birds
Meta <- fread("MetaData/Tagged bird summary data new.csv")
Meta <- filter(Meta, Ringing.location %in% c("ISLA", "WEST", "LOKE", "DYFI", "KINT") | S.N == 17795)

## filter the gps data for the device ids left
orn_gps <-filter(orn_gps, device_id %in% unique(Meta$S.N))

## parse timestamp
orn_gps$UTC_datetime <- as.character(ymd_hms(orn_gps$UTC_datetime))


## Delete observations where missing lat or long or a timestamp. 
ind <- complete.cases(orn_gps[, c("Latitude", "Longitude", "UTC_datetime")])
table(ind) # Number of re locations with missing coordinates or timestamp (if any)
orn_gps <- orn_gps %>% filter(ind) # removing rows with missing data (if any)

##Check for duplicated observations (ones with same lat, long, timestamp, and device ID)
ind2 <- orn_gps %>% 
        dplyr::select(UTC_datetime, Longitude, Latitude, device_id, acc_x, acc_z, acc_y) %>%
        duplicated()
sum(ind2) 
orn_gps <- orn_gps %>% filter(!ind2)


## Streamline the dataset for a first pass, first just have subset of the columns
orn_gps <- subset(orn_gps, select = c("device_id", "UTC_datetime", "Latitude", "Longitude"))
setnames(orn_gps, old = "device_id", new = "Tag_ID")

## Now filter for the winter
## create a year day column first
orn_gps$year_day <- yday(orn_gps$UTC_datetime)
winter_orn_gps <- filter(orn_gps, year_day >= 305 | year_day <= 92)
table(winter_orn_gps$Tag_ID)






##                                   ##
#### 1.3 Join together GPS data sets ####
##                                   ##

## bind the data sets
winter_eco_gps$UTC_datetime <- ymd_hms(winter_eco_gps$UTC_datetime)
winter_orn_gps$UTC_datetime <- ymd_hms(winter_orn_gps$UTC_datetime)
winter_gps <- rbind(winter_eco_gps, winter_orn_gps)


## Add a species column
winter_gps$species <- "GWfG"







##                                                  ##
#### 1.4 Read in GBG tracking and join to GWfG data ####
##                                                  ##

## read in the tracking data csv for GBG
GBG_GPS <- fread("tracking data/GBG_TAG_CLEAN.csv")

## extract columns that I want from the GBG data
GBG_sub <- subset(GBG_GPS, select = c("device_id", "UTC_datetime", "Longitude", "Latitude"))
## change column names of the GBG data so it matches the GWfG data
setnames(GBG_sub, old = c("device_id"), new = c("Tag_ID"))

## add a species column
GBG_sub$species <- "GBG"

## add year day column, parse timestamp first
GBG_sub$UTC_datetime <- as.POSIXct(GBG_sub$UTC_datetime, format = "%d/%m/%Y %H:%M:%S")
GBG_sub$year_day <- yday(ymd_hms(GBG_sub$UTC_datetime))

## filter out the winter period
GBG_winter <- filter(GBG_sub, year_day >= 305 | year_day <= 92)



## Now join to the GWfG data, make sure the timestamps are the same class
winter_gps$UTC_datetime <- ymd_hms(winter_gps$UTC_datetime)
GBG_sub$UTC_datetime <- ymd_hms(GBG_sub$UTC_datetime)
winter_gps  <- rbind(winter_gps, GBG_winter) 

##  order by TagID and by timestamp then add a unique ID to each fix
winter_gps <- winter_gps[order(winter_gps$Tag_ID, winter_gps$UTC_datetime),]
winter_gps$ID <- 1:nrow(winter_gps)




##
#### 2. Label with Landcover data ####
##

## All 5 years of land cover data have the same 21 habitat categories
## They make use of images throughout the FULL year to do the classification


##
#### 2.1 Read in raster data ####
##

## Read in the raster data 
LandCov_2015 <- raster("Landcover Data/Islay landcover data/LandCov2015_Islay.grd")
LandCov_2017 <- raster("Landcover Data/Islay landcover data/LandCov2017_Islay.grd")
LandCov_2018 <- raster("Landcover Data/Islay landcover data/LandCov2018_Islay.grd")
LandCov_2019 <- raster("Landcover Data/Islay landcover data/LandCov2019_Islay.grd")
LandCov_2020 <- raster("Landcover Data/Islay landcover data/LandCov2020_Islay.grd")

## Stack the Rasters
LandCov_stack <- raster::stack(LandCov_2015, LandCov_2017, LandCov_2018, LandCov_2019, LandCov_2020)
names(LandCov_stack) <- c("LC2015", "LC2017", "LC2018", "LC2019", "LC2020")




##
#### 2.2 Crop tracking data to same extent as that as Raster ####
##


## Convert the GPS data into a sf object
winter_gps_sf <- st_as_sf(winter_gps, coords = c("Longitude", "Latitude"), 
                          crs = 4326, agr = "constant")

## now convert to the same crs as the field data
winter_gps_sf <- st_transform(winter_gps_sf, crs = st_crs(LandCov_2015))

## crop to the same extent as the landcover data
Islay_gps_sf <- st_crop(winter_gps_sf, st_bbox(LandCov_2015))




##
#### 2.3 Add the ID of the land cover dataset each fix should be labeled with ####
##


Islay_gps_sf$LC <- ifelse(Islay_gps_sf$UTC_datetime < "2016-06-15 23:00:00", "LC2015", 
                          ifelse(Islay_gps_sf$UTC_datetime >= "2016-06-15 23:00:00" & Islay_gps_sf$UTC_datetime < "2018-06-15 23:00:00", "LC2017",
                                 ifelse(Islay_gps_sf$UTC_datetime >= "2018-06-15 23:00:00" & Islay_gps_sf$UTC_datetime < "2019-06-15 23:00:00", "LC2018", 
                                        ifelse(Islay_gps_sf$UTC_datetime >= "2019-06-15 23:00:00" & Islay_gps_sf$UTC_datetime < "2020-06-15 23:00:00", "LC2019", 
                                               ifelse(Islay_gps_sf$UTC_datetime >= "2020-06-15 23:00:00", "LC2020", "ERROR")))))

## check that row got assigned a land cover data set
table(Islay_gps_sf$LC)




##
#### 2.4 Label each fix with the correct landcover data set ####
##


## create a list of the unique values in the LC column
LandCovs <- unique(Islay_gps_sf$LC)

## now loop through each value in LandCovs and label with the corresponding landcover dataset
for(j in 1:length(LandCovs)){
  
  message(j, " out of ", length(LandCovs))
  
  ## extract the firsr round of fixes
  Sub_fixes <- filter(Islay_gps_sf, LC == LandCovs[j])
  
  ## extract the corresponding raster 
  Sub_raster <- subset(LandCov_stack, LandCovs[j])

  ## label with habitat from the corresponding raster
  Sub_fixes$Habitat <- raster::extract(x = Sub_raster, y = Sub_fixes, method = "simple")
  
  ## bind together all the sub datasets
  if(j ==1){All_fixes <- Sub_fixes}else{All_fixes <- rbind(All_fixes, Sub_fixes)}
  
}

## converting back to lat/longs
All_fixes_df <- st_transform(All_fixes, crs = 4326)
All_fixes_df <- sf_to_df(All_fixes_df, fill = T)


## remove columns added by sf
All_fixes_df$sfg_id <- NULL
All_fixes_df$point_id <- NULL
All_fixes_df$LC <- NULL
All_fixes_df$ID <- NULL

## set correct names for Lat/Long
setnames(All_fixes_df, old = c("x", "y"), new = c("Longitude", "Latitude"))


## save this data set to use in analysis
saveRDS(All_fixes_df, file = "Derived data/All_winter_GPS_with_habitat.RDS")







