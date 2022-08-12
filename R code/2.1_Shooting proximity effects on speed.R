## Luke Ozsanlav-Harris

## 13/01/2022

## Examining the effect that shooting has on movement, in particular the step length immediately after shooting disturbance. 
## Idea here is to model disturbed vs undisturbed step lengths
## This can either be done within individuals or between individuals
## Approach 1) within individuals. Compare step length before and during disturbance with an individual and compare
## Approach 2) between individual. Compare step lengths during disturbance to those from other individuals at the same time elsewere on the Island

# Approach 1) model: SL ~ disturbed + TOD + DOY^2 + (1|Pair ID) + (1|Tag ID) ..... these two random effects might need to be nested e.g. (1|Pair ID:Tag ID)
#                   note: might be easier to model the change between the control and treatment or a proportion of the two
# Approach 2) model:
#                   note: Might want to include the time difference between the control and treatment

## Packages required
library(amt)
library(tidyverse)
library(lubridate)
library(data.table)
library(zoo)
library(geosphere)
library(sf)
library(lme4)
library(DHARMa)
library(MuMIn)
library(effects)
library(performance)
library(lmerTest)
library(ggplot2)
library(MuMIn)

## Script Structure...
## 1.1 Read in and clean Ecotone GPS data
## 1.2 Read in and clean Ornitela GPS data
## 1.3 Join together GPS data sets
## 1.4 Read in GBG tracking and join to GWfG data
## 1.5 Read out GWfG tag codes for associations matrix
## 2.  Read in shooting data and field boundaries
## 3.  Combine shooting data and spatial field data
## 4.  Label fixes associated with each shooting event
## 5.  Deal with multi-associations and add shooting times
## 6.  Add shooting associations back to the full tracking data set
## 7.  Identify suitable paired observations within individuals
## 8.  Model how disturbance effects step length (Individual control)
## 9.  Plot how disturbance effects step length (Individual control)
## 10. Use threshold model to find distance were there is not shooting effect
## 11. Look at how the shooting response changes depending on which fix is chosen as the control


## Vignette for the effects package: https://cran.r-project.org/web/packages/effects/vignettes/predictor-effects-gallery.pdf

## Improvements overall:
## 2. Pick controls fixes from other birds that were undisturbed at the same time
## 3. plot output of effects package in ggplot
## 4. Can I use another package that will model species specific thresholds simultaneously and calculate if they are significantly different


## NOTES on results....
## Looking at species specific threshol responses
## From this it suggests that GWfG respond to shooting events further away but 
## that they have lower changed in speed siggesting they don't move as far as GBG in response to shooting

## I compared using the fix immediately before and after the disturbance fix. Didn't seem to influence the output at all....


##                                           ##
#### 1.1  Read in and clean Ecotone GPS data ####
##                                           ##


## set working directory
setwd("~/Ecotone Data/Cleaned Data file")

## read in data
## Will need to read in allthe data at a later point but just start with one file for now
#orn_data <-list.files(pattern = "*.csv") %>% map_df(~fread(.))
eco_data <- fread("Ecotone_all_clean.csv")


#### Data Cleaning 
## extracting just the GPS data
eco_gps <- filter(eco_data, is.na(Acc_Z) ==T)
rm(eco_data)

## parse timestamp
eco_gps$UTC_datetime <- as.character(ymd_hms(eco_gps$Date_Time_milli_posix))

## Delete observations where missing lat or long or a timestamp. 
colnames(eco_gps)
ind <- complete.cases(eco_gps[, c("Latitude", "Longitude", "UTC_datetime")])

## The number of re locations with missing coordinates or timestamp (if any).
table(ind)

##removing rows with missing data (if any)
eco_gps <- eco_gps %>% filter(ind)

##Check for duplicated observations (ones with same lat, long, timestamp, and device ID)
ind2 <- eco_gps %>% 
  dplyr::select(UTC_datetime, Longitude, Latitude, Tag_ID) %>%
  duplicated()
sum(ind2) 
eco_gps <- eco_gps %>% filter(!ind2)


## Streamline the dataset for a first pass
## first just have subset of the columns
colnames(eco_gps)
gps <- subset(eco_gps, select = c("Tag_ID" , "UTC_datetime", "Latitude", "Longitude"))

## Now filter for the winter
## create a year day column first
gps$year_day <- yday(gps$UTC_datetime)
winter_gps <- filter(gps, year_day >= 305 | year_day <= 92)
table(winter_gps$Tag_ID)







##                                            ##
#### 1.2  Read in and clean Ornitela GPS data ####
##                                            ##

## set working directory
setwd("~/Ornitela data/ALL_ORNI_DATA/Cleaned CSVs")

## read in data
## Will need to read in allthe data at a later point but just start with one file for now
#orn_data <-list.files(pattern = "*.csv") %>% map_df(~fread(.))
orn_data1 <- fread("SCOT1_clean.csv")
orn_data2 <- fread("SCOT2_clean.csv")
orndata <- rbind(orn_data1, orn_data2)
rm(orn_data1, orn_data2)

#### Data Cleaning
## extracting just the GPS data
orn_gps <- subset(orndata, orndata$datatype == "GPS")
rm(orndata)

## parse timestamp
orn_gps$UTC_datetime <- as.character(ymd_hms(orn_gps$UTC_datetime))


## Delete observations where missing lat or long or a timestamp. 
ind <- complete.cases(orn_gps[, c("Latitude", "Longitude", "UTC_datetime")])
table(ind) # The number of re locations with missing coordinates or timestamp (if any).
orn_gps <- orn_gps %>% filter(ind) # removing rows with missing data (if any)

##Check for duplicated observations (ones with same lat, long, timestamp, and device ID)
ind2 <- orn_gps %>% 
  dplyr::select(UTC_datetime, Longitude, Latitude, device_id, acc_x, acc_z, acc_y) %>%
  duplicated()
sum(ind2) 
orn_gps <- orn_gps %>% filter(!ind2)


## Streamline the dataset for a first pass
## first just have subset of the columns
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
winter_gps <- rbind(winter_gps, winter_orn_gps)
rm(winter_orn_gps)

## Add a species column
winter_gps$species <- "GWfG"







##                                                  ##
#### 1.4 Read in GBG tracking and join to GWfG data ####
##                                                  ##

## read in the tracking data csv for GBG
setwd("~/Shooting disturbance/tracking data")
GBG_GPS <- fread("GBG_TAG_CLEAN.csv")

## extract columns that i want from the GBG data
GBG_sub <- subset(GBG_GPS, select = c("device_id", "UTC_datetime", "Longitude", "Latitude"))
## change column names of the GBG data so it matches the GWfG data
setnames(GBG_sub, old = c("device_id"), new = c("Tag_ID"))

## add a species column
GBG_sub$species <- "GBG"

## add year day column
## parse timestamp first
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





##                                                       ##
#### 1.5 Read out GWfG tag codes for associations matrix ####
##                                                       ##

## Read in tag info so I can add sexes to the tags
setwd("~/Additional data files")
tag_info <- fread("Tagged bird summary data new.csv")

## Change serial numbers to character
tag_info$S.N <- as.character(tag_info$S.N)
## change the Bird IDs of the ornitela tags to the serial numbers
tag_info$Bird.ID <- ifelse(is.na(tag_info$S.N)==F, tag_info$S.N, tag_info$Bird.ID)

## take out just the columns I want
tag_info <- subset(tag_info, select = c("Bird.ID", "Age", "Sex"))

## Extract just the unique GWfG codes
GWfG_GPS <- filter(winter_gps, species == "GWfG")

## extract unique tag codes
codes <- as.data.frame(unique(GWfG_GPS$Tag_ID))
names(codes)[1] <- "Bird.ID"

## Add on tag info
codes2 <- left_join(codes, tag_info, by = "Bird.ID")

## write out the file
#write_csv(codes2, file = "~/Shooting disturbance/List of tags used/GWfG_winter_codes.csv")







##                                  ##
#### 2. Re sample all tracking data ####
##                                  ##


## Realized now (10/02/22) that i need to re-sample the tracking data to a standard sampling rate

## function to resample the tracking data
trackSubSampyear <- function(TD = TD, dt=dt, unit=unit){
  
  TD <- TD[order(TD$UTC_datetime),] 
  years = unique(TD$year)
  
  for(j in 1:length(years)){
    
    message(j, " out of ", length(years))
    
    TD_sub = filter(TD, year == years[j]) # extract one year at a time
    
    # breakdown to datasets per bird
    unid = unique(TD_sub$tag_year) 
    nrid = length(unid)
    TDall = list(nrid)  
    TDsubred = list(nrid)
    timestep = paste(dt,unit)
    # create time sequence from min to max time with step sizes that are defined at the start of the function
    dati_start = min(TD_sub$UTC_datetime, na.rm= T)
    dati_end = max(TD_sub$UTC_datetime, na.rm= T)
    datiseq = seq(from=dati_start, to=dati_end, by=timestep)
    
    for (i in 1:nrid) {
      Dtemp = TD_sub[TD_sub$tag_year == unid[i],]
      idx = sapply(datiseq, function(x) which.min( abs( difftime( Dtemp$UTC_datetime, x, units='mins')))) # finds closest time in data to your created time series
      TDall[[i]] = Dtemp
      TDsubred[[i]] = unique(Dtemp[idx,]) # the function unique makes sure that the rows in Dtemp[idx,] are unique - so no duplicate points
    }
    
    TDsubred2 <- do.call("rbind", TDsubred)
    if(j == 1){Resamp_GPS <- TDsubred2}
    else{Resamp_GPS <- rbind(Resamp_GPS, TDsubred2)}
    
  }
  
  return(Resamp_GPS)
  
}


## prepare data for function, timestamp needs to be POSIX object and add year column
winter_gps$UTC_datetime <- as.POSIXct(winter_gps$UTC_datetime, format= "%Y-%m-%d %H:%M:%OS", tz= "GMT")
winter_gps$year <- year(winter_gps$UTC_datetime)
winter_gps$tag_year <- paste0(winter_gps$Tag_ID, "_", winter_gps$year)

## re-sample the data, try 1 hour for now 
winter_resamp <- trackSubSampyear(TD = winter_gps, dt = 60, unit = 'mins')







##                                                 ##
#### 3. Read in shooting data and field boundaries ####
##                                                 ##

## read in the shooting logs
setwd("~/Shooting disturbance")
Logs <- fread("All_logs_cleaned.csv")

## create a timestamp column for the shooting logs 
## think the ones that fail to parse are the ones were the time is missing and it is just a date
Logs$timestamp <- ymd_hms(paste0(Logs$Date_Cl, " ", Logs$Time_Cl))


## Now filter the logs for the same time period as the shooting data
winter_logs <- filter(Logs, timestamp > min(winter_resamp$UTC_datetime) & timestamp < max(winter_resamp$UTC_datetime))

## Now filter out just the shooting events
winter_shots <- filter(winter_logs, Shots_fired_Cl > 0)

## read in the shapefiles of the Islay field boundaries
Fields <- st_read(dsn = "~/Shooting disturbance/Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY",
                  layer = "AS_ISLAY_GMS_FIELD_BOUNDARY")
st_crs(Fields)

## remove any duplicates in the field data
ind3 <- duplicated(Fields$FIELD_ID)
sum(ind3) 
Fields <- Fields %>% filter(!ind3)







##                                                   ##
#### 4. Combine shooting data and spatial field data ####
##                                                   ##


## calculate the centroid of each field
Field_centres <- st_centroid(Fields)

## now create a data set with the shooting events but where each one has the center for that field
setnames(winter_shots, old = "Field_code_Cl", new = "FIELD_ID")
winter_shots_cent <- left_join(winter_shots, Field_centres, by = "FIELD_ID")

## for now filter out the shooting events were we don't know the field code
winter_shots_cent2 <- filter(winter_shots_cent, is.na(SHAPE_LEN)==F)

## Give each shooting event a unique ID number
winter_shots_cent2$shot_ID <- 1:nrow(winter_shots_cent2)

## Make winter_shots_cent2 an sf object again
winter_shots_cent_sf <- st_as_sf(winter_shots_cent2)
winter_shots_cent_sf <- st_transform(winter_shots_cent_sf, crs = st_crs(Field_centres))







##                                                      ##
#### 5. Label fixes associated with each shooting event ####
##                                                      ##


## Define function that converts timestamp immediately before shooting event the smallest value
## y needs to be the timestamps from the tracking data and x is just a single value from the shooting data sequence
difftime2 <- function(x=x, y=y){
  
  time_difs <- difftime(y, x, units='mins')*-1
  time_difs[time_difs<0] <- 1000000000
  return(as.vector(time_difs))
}

## Define function that will loop through each tag and find the nearest fix before each shooting event while the tag was active
## The fix prior to the shooting event has to be within a certain spatial and temporal buffer to the shooting event
## TD is the tracking data. Need tag ID (Tag_ID), lat/long (Longitude/Latitude), timestamp (UTC_datetime) and unique ID for each for row (ID)
## SD is the shooting data. Need the time of the shot (timestamp), the field center (as an sf geometry) and unique ID for each for row (shot_ID)
## time_thresh is the max time in minutes that a fix can be before a shooting event
## dist_thresh is the max distance in meters the fix before a shooting event can be from the location of the shooting event
## field_crs is an sf object that has the crs of the field centers in SD

Shoot_proximity <- function(TD = TD, SD = SD, time_tresh = time_tresh, dist_thresh = dist_thresh, field_crs = field_crs){
  
  ##  make sure that the data set is ordered by time
  TD <- TD[order(TD$UTC_datetime),] 
  
  ## breakdown to datasets per bird
  unid <- unique(TD$Tag_ID) # list of unique tag IDs
  nrid <- length(unid) # number of unique birds
  
  ## Create lists to put each tag ID in
  Prox_list <- vector(mode = "list", length = nrid)
  
  ## Loop through each tag
  for (i in 1:nrid) { ## ** CHANGE 15 BACK TO nrid **###
    
    message("Tag ", i, " out of ", nrid)
    
    
    ## DATA PREP ##
    ## filter out each tag one at a time
    Dtemp <- TD[TD$Tag_ID == unid[i],]
    
    ## create sequence of shooting data while the tag was active
    Start_time <- min(Dtemp$UTC_datetime)
    End_time <- max(Dtemp$UTC_datetime)
    SD_sub <- filter(SD , timestamp > Start_time, timestamp < End_time)
    shootseq <- SD_sub$timestamp
    
    
    ## FINDING NEAREST PRIOR FIX ##
    # find the first fix before the shooting event
    idx = sapply(shootseq, function(x) which.min(difftime2(y = Dtemp$UTC_datetime, x = x)))
    
    ## combine these fixes immediately before the shooting event with the shooting data
    aa <- (Dtemp[idx,])
    bb <- cbind(aa , SD_sub)
    
    ##** COULD POTENTIALLY TAKE THIS SECTION OUT OF THE LOOP **##  
    
    ## CALCUALTE TIME DIFFERENCE AND DISTANCE ##
    ## calculate time differences between fix and shooting event
    bb$time_dif <- as.duration(ymd_hms(bb$timestamp) %--% ymd_hms(bb$UTC_datetime))
    bb$time_dif <- as.numeric(bb$time_dif, "minutes")
    
    ## Now calculate distances between the fix and the shooting event
    ## Make the fixes an sf object
    gps_sf <- st_as_sf(aa, coords = c("Longitude", "Latitude"), 
                       crs = 4326, agr = "constant")
    ## now convert to the same crs as the field data
    gps_sf <- st_transform(gps_sf, crs = st_crs(field_crs))
    
    ## make the shooting events an sf object
    shootseq_sf <- st_as_sf(SD_sub)
    shootseq_sf <- st_transform(shootseq_sf, crs = st_crs(field_crs))
    
    ## calculate the distance between the two sets of points
    distances <- st_distance(gps_sf, shootseq_sf, by_element = TRUE)
    
    ## put the distances back into the shooting proximity data set
    bb$dist_to_shot <- as.numeric(distances)
    
    
    ## RULES TO KEEP FIXES IN SPATIO TEMPORAL BUFFER ##
    ## Remove any fixes that were greater then 120 before the shooting event or were after the shooting event
    bb2 <- filter(bb, time_dif <= 0, time_dif >= -time_tresh)
    
    ## Remove any fixes that were greater than a certain number of meters from the shooting event
    bb3 <- filter(bb2, dist_to_shot < dist_thresh)
    
    ## Put this data set into a list that was initialized outside of the loop
    Prox_list[[i]] <- bb3
    
  }
  
  Prox_list
  
}



## run the function
## Used 2000m as the distance threshold, they ony tested up to c1,500m in the Islay study and at this distance the effect became small
## Added c500m on to allow a bit of wiggle room
system.time(Prox_list <- Shoot_proximity(TD = winter_resamp, SD = winter_shots_cent2, time_tresh = 120, dist_thresh = 4000, field_crs = Field_centres))

## bind the lists together into one data frame
All_prox <- plyr::ldply(Prox_list)








##                                                            ##
#### 6. Deal with multi associations and add shooting times ####
##                                                            ##

## some fixes overlap with multiple shooting events so need to sort that out first
## Will create a unique column for each shooting event

## Now combine all the lists into one data frame
All_prox <- plyr::ldply(Prox_list)
All_prox <- subset(All_prox, select = c("Tag_ID", "UTC_datetime", "ID", "Shots_fired_Cl", "timestamp", "shot_ID", "time_dif", "dist_to_shot"))

## Add a column that signifies these fixes overlap with shooting events
Shoot_prox <- All_prox
Shoot_prox$spatial_overlap <- 1

## summarize the data set for individuals fixes that spatio-temporally overlap with multiple shooting events
## First create function that extracts the shooting ID number but you can specify from which row of the data frame
ID_extractor <- function(x, row_num = row_num){
  
  return(as.numeric(x$shot_ID[row_num]))
  
}

time_extractor <- function(x, row_num = row_num){
  
  as.character(x$timestamp[row_num])
  
}

shots_extractor <- function(x, row_num = row_num){
  
  as.numeric(x$Shots_fired_Cl[row_num])
  
}

disttoshot_extractor <- function(x, row_num = row_num){
  
  as.numeric(x$dist_to_shot[row_num])
  
}

timedif_extractor <- function(x, row_num = row_num){
  
  as.numeric(x$time_dif[row_num])
  
}



## Now nest the data by GPS fix ID number
## For fixes associated with multiple shooting events there data frames will be > 1 row long
## Then just extract each shooting event in turn and put it in its own column
Multi_shoot <- Shoot_prox %>% 
  group_by(Tag_ID, UTC_datetime) %>% 
  nest() %>% 
  mutate(Shoot_ID1 = map(data, ID_extractor, row_num = 1),
         Shoot_ID2 = map(data, ID_extractor, row_num = 2), 
         Shoot_ID3 = map(data, ID_extractor, row_num = 3), 
         Shoot_ID4 = map(data, ID_extractor, row_num = 4), 
         Shoot_ID5 = map(data, ID_extractor, row_num = 5), 
         Shoot_ID1_time = map(data, time_extractor, row_num = 1),
         Shoot_ID2_time = map(data, time_extractor, row_num = 2), 
         Shoot_ID3_time = map(data, time_extractor, row_num = 3), 
         Shoot_ID4_time = map(data, time_extractor, row_num = 4), 
         Shoot_ID5_time = map(data, time_extractor, row_num = 5),
         Shoot_ID1_shots = map(data, shots_extractor, row_num = 1),
         Shoot_ID2_shots = map(data, shots_extractor, row_num = 2), 
         Shoot_ID3_shots = map(data, shots_extractor, row_num = 3), 
         Shoot_ID4_shots = map(data, shots_extractor, row_num = 4), 
         Shoot_ID5_shots = map(data, shots_extractor, row_num = 5),
         Shoot_ID1_dist_to = map(data, disttoshot_extractor, row_num = 1),
         Shoot_ID2_dist_to = map(data, disttoshot_extractor, row_num = 2), 
         Shoot_ID3_dist_to = map(data, disttoshot_extractor, row_num = 3), 
         Shoot_ID4_dist_to = map(data, disttoshot_extractor, row_num = 4), 
         Shoot_ID5_dist_to = map(data, disttoshot_extractor, row_num = 5),
         Shoot_ID1_time_dif = map(data, timedif_extractor, row_num = 1),
         Shoot_ID2_time_dif = map(data, timedif_extractor, row_num = 2), 
         Shoot_ID3_time_dif = map(data, timedif_extractor, row_num = 3), 
         Shoot_ID4_time_dif = map(data, timedif_extractor, row_num = 4), 
         Shoot_ID5_time_dif = map(data, timedif_extractor, row_num = 5))

## remove any columns that only contain NAs
## first drop the data column, then remove columns with only NAs
Multi_shoot <- Multi_shoot %>% 
  dplyr::select(c(-"data")) %>% 
  as.data.frame() %>% 
  select_if(~sum(!is.na(.)) > 0)

## set the shooting ID columns to numeric and time to lubridate objects
Multi_shoot$Shoot_ID1 <- as.numeric(Multi_shoot$Shoot_ID1)
Multi_shoot$Shoot_ID2 <- as.numeric(Multi_shoot$Shoot_ID2)
Multi_shoot$Shoot_ID3 <- as.numeric(Multi_shoot$Shoot_ID3)
Multi_shoot$Shoot_ID4 <- as.numeric(Multi_shoot$Shoot_ID4)
Multi_shoot$Shoot_ID5 <- as.numeric(Multi_shoot$Shoot_ID5)
Multi_shoot$Shoot_ID1_time <- ymd_hms(Multi_shoot$Shoot_ID1_time)
Multi_shoot$Shoot_ID2_time <- ymd_hms(Multi_shoot$Shoot_ID2_time)
Multi_shoot$Shoot_ID3_time <- ymd_hms(Multi_shoot$Shoot_ID3_time)
Multi_shoot$Shoot_ID4_time <- ymd_hms(Multi_shoot$Shoot_ID4_time)
Multi_shoot$Shoot_ID5_time <- ymd_hms(Multi_shoot$Shoot_ID5_time)
Multi_shoot$Shoot_ID1_shots <- as.numeric(Multi_shoot$Shoot_ID1_shots)
Multi_shoot$Shoot_ID2_shots <- as.numeric(Multi_shoot$Shoot_ID2_shots)
Multi_shoot$Shoot_ID3_shots <- as.numeric(Multi_shoot$Shoot_ID3_shots)
Multi_shoot$Shoot_ID4_shots <- as.numeric(Multi_shoot$Shoot_ID4_shots)
Multi_shoot$Shoot_ID5_shots <- as.numeric(Multi_shoot$Shoot_ID5_shots)
Multi_shoot$Shoot_ID1_dist_to <- as.numeric(Multi_shoot$Shoot_ID1_dist_to)
Multi_shoot$Shoot_ID2_dist_to <- as.numeric(Multi_shoot$Shoot_ID2_dist_to)
Multi_shoot$Shoot_ID3_dist_to <- as.numeric(Multi_shoot$Shoot_ID3_dist_to)
Multi_shoot$Shoot_ID4_dist_to <- as.numeric(Multi_shoot$Shoot_ID4_dist_to)
Multi_shoot$Shoot_ID5_dist_to <- as.numeric(Multi_shoot$Shoot_ID5_dist_to)
Multi_shoot$Shoot_ID1_time_dif <- as.numeric(Multi_shoot$Shoot_ID1_time_dif)
Multi_shoot$Shoot_ID2_time_dif <- as.numeric(Multi_shoot$Shoot_ID2_time_dif)
Multi_shoot$Shoot_ID3_time_dif <- as.numeric(Multi_shoot$Shoot_ID3_time_dif)
Multi_shoot$Shoot_ID4_time_dif <- as.numeric(Multi_shoot$Shoot_ID4_time_dif)
Multi_shoot$Shoot_ID5_time_dif <- as.numeric(Multi_shoot$Shoot_ID5_time_dif)


## Remove duplicate rows from the main data set
ind_shot <- duplicated(Shoot_prox$ID)
table(ind_shot)
Shoot_prox <- filter(Shoot_prox, !ind_shot)

## Now add on the multiple columns created for multiple shooting events
Shoot_prox <- left_join(Shoot_prox, Multi_shoot, by = c("Tag_ID", "UTC_datetime"))
Shoot_prox$shot_ID <- NULL # no longer need this column
Shoot_prox$timestamp <- NULL
Shoot_prox$Shots_fired_Cl <- NULL
Shoot_prox$dist_to_shot <- NULL
Shoot_prox$time_dif <- NULL








##                                                                   ##
#### 7. Add shooting associations back to the full tracking data set ####
##                                                                   ##


## Join the proximity data set back to the full gps set
colnames(Shoot_prox)
Shoot_prox2 <- subset(as.data.frame(Shoot_prox), select = -c(Tag_ID, UTC_datetime))
winter_gps_prox <- full_join(winter_resamp, Shoot_prox2, by = "ID")

## assign all NAs in the spatial overlap column a 0
winter_gps_prox$spatial_overlap <- ifelse(is.na(winter_gps_prox$spatial_overlap) == T, 0, winter_gps_prox$spatial_overlap)


#### Calculate step lengths and time between the fixes ####

## create a distance calculator function
## calculates distance moved since the last fix
Distance_calc <- function(x){
  
  distGeo(cbind(x$Longitude, x$Latitude), 
          cbind(lead(x$Longitude), lead(x$Latitude)))/1000
  
}

Time_dif_cal <- function(x){
  
  #as.numeric(lead(ymd_hms(x$UTC_datetime)) - ymd_hms(x$UTC_datetime))
  TD <- as.duration(ymd_hms(x$UTC_datetime) %--% lead(ymd_hms(x$UTC_datetime)))
  as.numeric(TD, "minutes")
  
}

## create a tag year column
winter_gps_prox$tag_year <- paste0(winter_gps_prox$Tag_ID, "_", year(winter_gps_prox$UTC_datetime))

## now nest the data frame to calculate the step lengths within tag years
winter_gps_prox <- winter_gps_prox %>% 
                    group_by(tag_year) %>% 
                    nest() %>% 
                    mutate(Step_len = map(data, Distance_calc),
                           time_dif = map(data, Time_dif_cal)) %>% 
                    unnest(cols = c(data, Step_len, time_dif))







##                                                               ##
#### 8. Identify suitable paired observations within individuals ####
##                                                               ##

##
#### 8.1 Shoot and t-1 as control ####
##

## Compare the disturbed fix to both the fix immediately before and after
## Copy the main data set just to find suitable pairs for now
pairs <- winter_gps_prox


## Add columns to the dataset for the time difs, shoot status and ID of the previous row....
ppairs <- pairs
ppairs$TD_control <- lag(ppairs$time_dif)
ppairs$shoot_control <- ifelse(is.na(lag(ppairs$Shoot_ID1))==T, 0, 1)
ppairs$ID_control <- lag(ppairs$ID)
ppairs$SL_control <- lag(ppairs$Step_len)

## Now just filter out the rows that experienced a shooting event, to streamline data set a bit now
ppairs2 <- ppairs %>% filter(is.na(Shoot_ID1) == F) 
## filter out instances were the disturbance step length is calculated from a fix interval greater than 130 mins
ppairs2 <- ppairs2 %>% filter(time_dif < 130)
## filter out instances when the fix interval of the proposed pair is greater than 10 minutes different
ppairs2 <- ppairs2 %>% filter(TD_control < time_dif + 10 & TD_control > time_dif - 10)
## filter out instances where the paired observation is also a disturbed
ppairs2 <- ppairs2 %>% filter(shoot_control == 0)



## Now extract the pairs and give each pair and their time dif and step length
Dist_fix <- subset(ppairs2, select = c("ID", "time_dif", "Step_len"))
Pair_fix <- subset(ppairs2, select = c("ID_control", "TD_control", "SL_control"))

## label each pair
Dist_fix$Pair_ID <- 1:nrow(ppairs2)

## calculate the speeds for both data sets, swapped these around as the 
Dist_fix$speed <- (Dist_fix$Step_len*1000)/(Dist_fix$time_dif*60) #meters/seconds
Pair_fix$prev_speed <- (Pair_fix$SL_control*1000)/(Pair_fix$TD_control*60)

## bind the data sets side by side now
final_pairs <- cbind(Dist_fix, Pair_fix)

## signify whether it is before of after the disturbance fix
final_pairs$control <- "Shoot_T-1"




##
#### 8.2 Shoot and t+1 as control ####
##

## Join on the time difs, shoot status and ID of the next row....
npairs <- pairs
npairs$TD_control <- lead(npairs$time_dif)
npairs$shoot_control <- ifelse(is.na(lead(npairs$Shoot_ID1))==T, 0, 1)
npairs$ID_control <- lead(npairs$ID)
npairs$SL_control <- lead(npairs$Step_len)


## Now just filter out the rows that experienced a shooting event, to streamline data set a bit now
npairs2 <- npairs %>% filter(is.na(Shoot_ID1) == F) 
## filter out instances were the disturbance step length is calculated from a fix interval greater than 130 mins
npairs2 <- npairs2 %>% filter(time_dif < 130)
## filter out instances when the fix interval of the proposed pair is greater than 10 minutes different
npairs2 <- npairs2 %>% filter(TD_control < time_dif + 10 & TD_control > time_dif - 10)
## filter out instances where the paired observation is also a disturbed
npairs2 <- npairs2 %>% filter(shoot_control == 0)



## Now extract the pairs and give each pair and their time dif and step length
Dist_fix2 <- subset(npairs2, select = c("ID", "time_dif", "Step_len"))
Pair_fix2 <- subset(npairs2, select = c("ID_control", "TD_control", "SL_control"))

## label each pair
Dist_fix2$Pair_ID <- 1:nrow(npairs2)

## calculate the speeds for both data sets
Dist_fix2$speed <- (Dist_fix2$Step_len*1000)/(Dist_fix2$time_dif*60) #meters/seconds
Pair_fix2$prev_speed <- (Pair_fix2$SL_control*1000)/(Pair_fix2$TD_control*60)

## bind the data sets side by side now
final_pairs2 <- cbind(Dist_fix2, Pair_fix2)

## signify whether it is before of after the disturbance fix
final_pairs2$control <- "Shoot_T+1"



##
#### 8.3 t+1 and t-1 as control ####
##

## Join on the time difs, shoot status and ID of the next row....
## want t-1 to be the control and t+1 to be the comparison
cpairs <- pairs
cpairs$TD_control <- lag(cpairs$time_dif, n = 2L)
cpairs$shoot_control <- ifelse(is.na(lag(cpairs$Shoot_ID1, n = 2L))==T, 0, 1)
cpairs$ID_control <- lag(cpairs$ID, n = 2L)
cpairs$SL_control <- lag(cpairs$Step_len, n = 2L)
cpairs$Tminus1 <- ifelse(is.na(lag(cpairs$Shoot_ID1))==T, 0, 1)
cpairs$ID_Shot <- ifelse(is.na(lag(cpairs$Shoot_ID1))==F, lag(cpairs$ID), NA) # ID of the shooting event 


## Now just filter out the rows that are immediately prior to a shooting event (t-1), to streamline data set a bit now
cpairs2 <- cpairs %>% filter(Tminus1 == 1) 
##filter out the data if the t-1 step length was also exposed to shooting
cpairs2 <- cpairs2 %>% filter(is.na(Shoot_ID1) == T)
## filter out instances were the (t-1) step length is calculated from a fix interval greater than 130 mins
cpairs2 <- cpairs2 %>% filter(time_dif < 130)
## filter out instances when the fix interval of the proposed pair is greater than 10 minutes different
cpairs2 <- cpairs2 %>% filter(TD_control < time_dif + 10 & TD_control > time_dif - 10)
## filter out instances where the paired observation is also a disturbed
cpairs2 <- cpairs2 %>% filter(shoot_control == 0)



## Now extract the pairs and give each pair and their time dif and step length
Dist_fix3 <- subset(cpairs2, select = c("ID", "time_dif", "Step_len", "ID_Shot"))
Pair_fix3 <- subset(cpairs2, select = c("ID_control", "TD_control", "SL_control"))

## label each pair
Dist_fix3$Pair_ID <- 1:nrow(cpairs2)

## calculate the speeds for both data sets
Dist_fix3$speed <- (Dist_fix3$Step_len*1000)/(Dist_fix3$time_dif*60) #meters/seconds
Pair_fix3$prev_speed <- (Pair_fix3$SL_control*1000)/(Pair_fix3$TD_control*60)

## bind the data sets side by side now
final_pairsX <- cbind(Dist_fix3, Pair_fix3)

## signify whether it is before of after the disturbance fix
final_pairsX$control <- "T-1_T+1"



##
#### 8.4 Bind all pairs together ####
##

## Bind the before and after paris together
final_pairs3 <- plyr::rbind.fill(final_pairs, final_pairs2, final_pairsX)

## for t-1/t+1 we had to already assign the ID if the shooting event fix as it isnt either of the fixes in the sample
## for the other tow samples the shooting event ID is just the ID
final_pairs3$ID_Shot <- ifelse(is.na(final_pairs3$ID_Shot)==T, final_pairs3$ID, final_pairs3$ID_Shot)

## remove the columns that are already present in final_pairs3
winter_gps_prox2 <- subset(winter_gps_prox, select = c(-Step_len, -time_dif))

## bind to final_pairs data set
final_pairs4 <- left_join(final_pairs3, winter_gps_prox2, by = c("ID_Shot" = "ID"))
stopifnot(nrow(final_pairs4) == nrow(final_pairs3))

## Just going to run this analysis with pairs were the fix length is roughly one hour
final_pairs4 <- filter(final_pairs4, time_dif > 50 & time_dif < 70)





##                                                                   
#### 9. Prepare data for models ####
##                                                                    

## Some disturbance fixes are exposed to multiple shooting events
## Therefore going to pick the mimimum distance for each for
## First pick out the column with the distanc eto shot data
Shot_dists <- subset(final_pairs4, select =c("Shoot_ID1_dist_to", "Shoot_ID2_dist_to", "Shoot_ID3_dist_to", "Shoot_ID4_dist_to"))
final_pairs4$min_shot_dist <- apply(Shot_dists, 1, FUN = min, na.rm = TRUE)
summary(final_pairs4$min_shot_dist)

## If I modeled these as differences in speeds then this would standardize for fix length
final_pairs4$speed_dif <- final_pairs4$speed - final_pairs4$prev_speed


## Add other variables to the data set that will go in the model
## year and column
final_pairs4$year <- year(final_pairs4$UTC_datetime)
## shooting time of day column
final_pairs4$TOD_ <- as.numeric(hour(final_pairs4$Shoot_ID1_time)) + as.numeric(minute(final_pairs4$Shoot_ID1_time)/60)
## days Since November 1st column
final_pairs4$year_day <- yday(final_pairs4$UTC_datetime)
final_pairs4$from_nov1 <- ifelse(final_pairs4$year_day > 250, 
                           final_pairs4$year_day - 305, 
                           final_pairs4$year_day + 60)
## add winter column
final_pairs4$winter <- ifelse(final_pairs4$year_day < 150, paste0((final_pairs4$year-1), "-", final_pairs4$year), paste0(final_pairs4$year, "-", (final_pairs4$year+1)))


## scale continuous variables
final_pairs4$TOD_sc <- scale(final_pairs4$TOD_)
final_pairs4$from_nov1_sc <- scale(final_pairs4$from_nov1)
final_pairs4$min_shot_dist_sc <- scale(final_pairs4$min_shot_dist)

## set variables as the right classes
final_pairs4$Tag_ID <- as.factor(final_pairs4$Tag_ID)
final_pairs4$year <- as.factor(final_pairs4$year)
final_pairs4$species <- as.factor(final_pairs4$species)

## check the distribution of the response
hist(final_pairs4$speed_dif)


  





##                                                                              
#### 10. Use threshold model to find distance were there is not shooting effect ####
##  

##
#### 10.1 Define functions for threshold model ####
##

## Single Threshold Model Function
## This goes into the lmer model formula
threshold1<-function(age,T1)
{
  # T1 threshold for age
  #
  age.1 <- age.2 <- age
  age.1[age.1 > T1] <- T1
  age.2[age.2 <= T1] <- 0
  age.2[age.2 > T1] <- age.2[age.2 > T1] - T1
  cbind(age1.1 = age.1, age1.2 = age.2)
}


# Make empty outputs to store model AICs
summary(final_pairs4$min_shot_dist)
output <- data.frame(model=as.numeric(), AIC=numeric()) # This is for single threshold model 

## Pick which data set I want to put inot the models
unique(final_pairs4$control)
Consec_pair <- filter(final_pairs4, control == "Shoot_T-1")
Consec_pair <- filter(Consec_pair, is.na(species) == F) # removes any rows were the data for the shooting event is missing
Gap_pair <- filter(final_pairs4, control == "T-1_T+1")
Gap_pair <- filter(Gap_pair, is.na(species) == F) # removes any rows were the data for the shooting event is missing



##
#### 11.1: Option 1: model both species together, forcing them to have same threshold ####
##

## create sequence of distances, each of which will be used as the break point in the model
dist_seq <- seq(0,3950,50)

## create its own output data set 
output1 <- output

## loop through each threshold value individual and save model AICc in output data set
## List of possible predictors: species + TOD_sc + winter + from_nov1_sc + I(from_nov1_sc^2) + (1|Shoot_ID1) + (1|Tag_ID)

## start at 2 as 1 is 0metres
for(i in 2:length(dist_seq)) {
  
  message(i, " out of ", length(dist_seq))
  
  model <- lmer(speed_dif ~ threshold1(min_shot_dist, dist_seq[i])*species + winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                data=Consec_pair,
                REML=FALSE)
  
  output1[i,]<-c(dist_seq[i], MuMIn::AICc(model))
  
}

output1
plot(output1$model, output1$AIC)
dist_seq[which.min(output1$AIC)] ## threshold point in metres


#### Plot output of the top threshold model
top.model <- lmer(speed_dif ~ threshold1(min_shot_dist, 1350)*species + winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                  data=Consec_pair,
                  REML=FALSE)

## use the effects package, and ggplot to plot the model predictions
## first the effects for each predicitor
divisions <- 200
top_mod_effects <- predictorEffects(top.model, focal.levels = divisions)
model_performance(top.model)

## now extract the fits for the first variable and bind them together
effects1 <- top_mod_effects[1]
fit1 <- as.data.frame(cbind(effects1[["min_shot_dist"]][["fit"]], effects1[["min_shot_dist"]][["lower"]], 
                            effects1[["min_shot_dist"]][["upper"]], effects1[["min_shot_dist"]][["x"]][["min_shot_dist"]]))
## change the names to something meaningful
setnames(fit1, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "dist_to_shot"))
## add a species column first half is GBG as that come first in the alphabet
fit1$Species <- NA
fit1$Species[1:divisions] <- "Barnacle"
fit1$Species[(divisions+1):(divisions*2)] <- "White-fronted"

## Now plot using ggplot
ggplot() + 
  geom_line(data=fit1, aes(x= dist_to_shot, y = fit, group = Species, colour = Species), size = 1.25)  +
  geom_ribbon(data = fit1, aes(x=dist_to_shot, ymin = lower, ymax = upper, group = Species, colour = Species), 
              alpha = 0.3, colour = NA, fill = "grey") + 
  ylab("Speed difference (m/s)") + xlab("Distance to shot (m)") +
  theme_bw() +
  scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  xlim (0, 3000) +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank())








##
#### 11.2 Option 2: model both species separately, Different thresholds ####
##

## create its own output data set 
outputGBG <- output
outputGWfG <- output

## filter the data set for each species
GBG_pairs <- filter(Consec_pair, species == "GBG")
GWfG_pairs <- filter(Consec_pair, species == "GWfG")


##
#### GBG threshold model
##

## loop through each threshold value individual and save model AICc in ouput data set
for(i in 2:length(dist_seq)) {
  
  message(i, " out of ", length(dist_seq))
  
  model <- lmer(speed_dif ~ threshold1(min_shot_dist, dist_seq[i]) + winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                data= GBG_pairs,
                REML=FALSE)
  
  outputGBG[i,]<-c(dist_seq[i], MuMIn::AICc(model))
  
}

outputGBG
ggplot(data = outputGBG) + geom_point(aes(x= model, y= AIC)) +
  xlab("Distance to shot (m)") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank())
setwd("~/Shooting disturbance/Plots for thesis chapter/Script 2) plots")
ggsave("AIC GBG single threshold.png", 
       width = 20, height = 15, units = "cm")
dist_seq[which.min(outputGBG$AIC)] # threshold point in meters


#### Plot output of the top threshold model
top.modelGBG <- lmerTest::lmer(speed_dif ~ threshold1(min_shot_dist, 1000) + winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                  data=GBG_pairs,
                  REML=FALSE)
summary(top.modelGBG)



## use the effects package to extract the fit for the first variable
GBG_mod_effects <- predictorEffects(top.modelGBG, focal.levels = divisions)
plot(GBG_mod_effects[1])
effectsGBG <- GBG_mod_effects[1]
fitGBG <- as.data.frame(cbind(effectsGBG[["min_shot_dist"]][["fit"]], effectsGBG[["min_shot_dist"]][["lower"]], 
                            effectsGBG[["min_shot_dist"]][["upper"]], effectsGBG[["min_shot_dist"]][["x"]][["min_shot_dist"]]))
## change the names to something meaningful
setnames(fitGBG, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "dist_to_shot"))
## add a species column first half is GBG as that come first in the alphabet
fitGBG$Species <- "Barnacle"



##
#### GWfG threshold model
##

## loop through each threshold value individual and save model AICc in ouput data set
for(i in 2:length(dist_seq)) {
  
  message(i, " out of ", length(dist_seq))
  
  model <- lmer(speed_dif ~ threshold1(min_shot_dist, dist_seq[i]) +  winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                data= GWfG_pairs,
                REML=FALSE)
  
  outputGWfG[i,]<-c(dist_seq[i], MuMIn::AICc(model))
  
}

outputGWfG
ggplot(data = outputGWfG) + geom_point(aes(x= model, y= AIC)) +
  xlab("Distance to shot (m)") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank())
setwd("~/Shooting disturbance/Plots for thesis chapter/Script 2) plots")
ggsave("AIC GWfG single threshold.png", 
       width = 20, height = 15, units = "cm")
dist_seq[which.min(outputGWfG$AIC)] # threshold point in meters


#### Plot output of the top threshold model
## use lmerTest::lmer to get p-values using the Satterthwaite approximation
top.modelGWfG <- lmerTest::lmer(speed_dif ~ threshold1(min_shot_dist, 1350) + winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                     data=GWfG_pairs,
                     REML=FALSE)
summary(top.modelGWfG)
drop1(top.modelGWfG)
confint(top.modelGWfG)
MuMIn::r.squaredGLMM(top.modelGWfG)


## use the effects package to extract the fit for the first variable
GWfG_mod_effects <- predictorEffects(top.modelGWfG, focal.levels = divisions)
plot(GWfG_mod_effects[1])
effectsGWfG <- GWfG_mod_effects[1]
fitGWfG <- as.data.frame(cbind(effectsGWfG[["min_shot_dist"]][["fit"]], effectsGWfG[["min_shot_dist"]][["lower"]], 
                              effectsGWfG[["min_shot_dist"]][["upper"]], effectsGWfG[["min_shot_dist"]][["x"]][["min_shot_dist"]]))
## change the names to something meaningful
setnames(fitGWfG, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "dist_to_shot"))
## add a species column first half is GWfG as that come first in the alphabet
fitGWfG$Species <- "White-fronted"



#### Plot the output of the separate GBG and GWfG models on the same plot
fit2 <- rbind(fitGBG, fitGWfG)
## Now plot using ggplot
ggplot() + 
  geom_line(data=fit2, aes(x= dist_to_shot, y = fit, group = Species, colour = Species), size = 1.25)  +
  geom_ribbon(data = fit2, aes(x=dist_to_shot, ymin = lower, ymax = upper, group = Species, colour = Species), 
              alpha = 0.3, colour = NA, fill = "grey") + 
  ylab("Speed difference (m/s)") + xlab("Distance to shot (m)") +
  theme_bw() +
  scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  xlim (0, 3000) +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank())






##
#### 11.3 Option 3: model GBGs separately, two thresholds points ####
##

## create thresholds to trial
dist_seq2 <- expand.grid(t1 = seq(50,3950,50), t2 =seq(50,3950,50))
dist_seq2 <- filter(dist_seq2, t1 < t2)

## create its own output data set 
outputGBG2 <- output

## filter the data set for each species
GBG_pairs2 <- filter(Consec_pair, species == "GBG")


##
#### GBG double threshold model
##

## loop through each threshold value individual and save model AICc in ouput data set
for(i in 2:nrow(dist_seq2)) {
  
  svMisc::progress(i, max.value = nrow(dist_seq2))
  
  model <- lmer(speed_dif ~ threshold2(min_shot_dist, dist_seq2$t1[i], dist_seq2$t2[i]) + winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                data= GBG_pairs2,
                REML=FALSE)
  
  outputGBG2[i,]<-c(paste0(dist_seq2$t1[i], " & ", dist_seq2$t2[i]), MuMIn::AICc(model))
  
}

outputGBG2
plot(outputGBG2$model, outputGBG2$AIC)
dist_seq2[which.min(outputGBG2$AIC),] # threshold point in meters


#### Plot output of the top threshold model
top.modelGBG2 <- lmerTest::lmer(speed_dif ~ threshold2(min_shot_dist, 150, 1050) + winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                               data=GBG_pairs2,
                               REML=FALSE)
summary(top.modelGBG2)
drop1(top.modelGBG2)
confint(top.modelGBG2)
MuMIn::r.squaredGLMM(top.modelGBG2)


## use the effects package to extract the fit for the first variable
GBG_mod_effects2 <- predictorEffects(top.modelGBG2, focal.levels = divisions)
plot(GBG_mod_effects2[1])
effectsGBG2 <- GBG_mod_effects2[1]
fitGBG2 <- as.data.frame(cbind(effectsGBG2[["min_shot_dist"]][["fit"]], effectsGBG2[["min_shot_dist"]][["lower"]], 
                              effectsGBG2[["min_shot_dist"]][["upper"]], effectsGBG2[["min_shot_dist"]][["x"]][["min_shot_dist"]]))
## change the names to something meaningful
setnames(fitGBG2, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "dist_to_shot"))
## add a species column first half is GBG as that come first in the alphabet
fitGBG2$Species <- "Barnacle"


ggplot() + 
  geom_line(data=fitGBG2, aes(x= dist_to_shot, y = fit, colour = Species), size = 1.25)  +
  geom_ribbon(data = fitGBG2, aes(x=dist_to_shot, ymin = lower, ymax = upper), 
              alpha = 0.3, colour = NA, fill = "grey") + 
  ylab("Speed difference (m/s)") + xlab("Distance to shot (m)") +
  theme_bw() +
  scale_colour_manual(values=c("#0072B2")) +
  xlim (0, 3000) +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank())


## Find the point the lower CI equals zero, going to use this as out cut off shor disturbed/undisturbed
top_mod_effects <- predictorEffects(top.modelGBG2, focal.levels = 10000)
lowerCI <- top_mod_effects[["min_shot_dist"]][["lower"]] ## lower CI estimates
Distances <- top_mod_effects[["min_shot_dist"]][["x"]][["min_shot_dist"]] ## distance from shot series
Distances[which.min(abs(lowerCI))] # distance were lower CI is closest to zero


##
#### 11.4 Plot the best models for t and t-1 ####
##

## The best models were the single threshold model for GWfG and double threshold for GBG
Bestfits <- rbind(fitGBG2, fitGWfG) ##vreate dats set

Bestfits$Species <- ifelse(Bestfits$Species == "White-fronted", "GWfG", "GBG")


ggplot() + 
  geom_line(data=Bestfits, aes(x= dist_to_shot, y = fit, group = Species, colour = Species), size = 1.25)  +
  geom_ribbon(data = Bestfits, aes(x=dist_to_shot, ymin = lower, ymax = upper, group = Species, colour = Species), 
              alpha = 0.3, colour = NA, fill = "grey") + 
  ylab("Speed difference (m/s) [t-(t-1)]") + xlab("Distance to shot (m)") +
  theme_bw() +
  scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  xlim (0, 3000) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title=element_text(size=12,), 
        axis.text = element_text(size =14), 
        axis.title.x = element_text(size =16),
        axis.title.y = element_text(size =16), 
        strip.text.x = element_text(size =14), 
        legend.title = element_text(size =14),
        legend.text = element_text(size =14))






setwd("~/Shooting disturbance/Plots for thesis chapter/Script 2) plots")
ggsave("Displacement as function of distance to shot [t-(t-1)].png", 
        width = 25, height = 15, units = "cm")










##
#### 11.5 Option 4: use t-1 and t+1 ####
##

## create its own output data set 
outputGBG <- output
outputGWfG <- output

## filter the data set for each species
GBG_pairs <- filter(Gap_pair, species == "GBG")
GWfG_pairs <- filter(Gap_pair, species == "GWfG")


##
#### GBG threshold model
##

## loop through each threshold value individual and save model AICc in ouput data set
for(i in 2:length(dist_seq)) {
  
  message(i, " out of ", length(dist_seq))
  
  model <- lmer(speed_dif ~ threshold1(min_shot_dist, dist_seq[i]) + winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                data= GBG_pairs,
                REML=FALSE)
  
  outputGBG[i,]<-c(dist_seq[i], MuMIn::AICc(model))
  
}

outputGBG
plot(outputGBG$model, outputGBG$AIC)
dist_seq[which.min(outputGBG$AIC)] # threshold point in meters


#### Plot output of the top threshold model
top.modelGBG <- lmerTest::lmer(speed_dif ~ threshold1(min_shot_dist, 400) + winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                               data=GBG_pairs,
                               REML=FALSE)
summary(top.modelGBG)


## use the effects package to extract the fit for the first variable
GBG_mod_effects <- predictorEffects(top.modelGBG, focal.levels = divisions)
plot(GBG_mod_effects[1])
effectsGBG <- GBG_mod_effects[1]
fitGBG <- as.data.frame(cbind(effectsGBG[["min_shot_dist"]][["fit"]], effectsGBG[["min_shot_dist"]][["lower"]], 
                              effectsGBG[["min_shot_dist"]][["upper"]], effectsGBG[["min_shot_dist"]][["x"]][["min_shot_dist"]]))
## change the names to something meaningful
setnames(fitGBG, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "dist_to_shot"))
## add a species column first half is GBG as that come first in the alphabet
fitGBG$Species <- "Barnacle"



##
#### GWfG threshold model
##

## loop through each threshold value individual and save model AICc in ouput data set
for(i in 2:length(dist_seq)) {
  
  message(i, " out of ", length(dist_seq))
  
  model <- lmer(speed_dif ~ threshold1(min_shot_dist, dist_seq[i]) +  winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                data= GWfG_pairs,
                REML=FALSE)
  
  outputGWfG[i,]<-c(dist_seq[i], MuMIn::AICc(model))
  
}

outputGWfG
plot(outputGWfG$model, outputGWfG$AIC)
dist_seq[which.min(outputGWfG$AIC)] # threshold point in meters


#### Plot output of the top threshold model
## use lmerTest::lmer to get p-values using the Satterthwaite approximation
top.modelGWfG <- lmerTest::lmer(speed_dif ~ threshold1(min_shot_dist, 1350) + winter + TOD_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                                data=GWfG_pairs,
                                REML=FALSE)
summary(top.modelGWfG)


## use the effects package to extract the fit for the first variable
GWfG_mod_effects <- predictorEffects(top.modelGWfG, focal.levels = divisions)
plot(GWfG_mod_effects[1])
effectsGWfG <- GWfG_mod_effects[1]
fitGWfG <- as.data.frame(cbind(effectsGWfG[["min_shot_dist"]][["fit"]], effectsGWfG[["min_shot_dist"]][["lower"]], 
                               effectsGWfG[["min_shot_dist"]][["upper"]], effectsGWfG[["min_shot_dist"]][["x"]][["min_shot_dist"]]))
## change the names to something meaningful
setnames(fitGWfG, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "dist_to_shot"))
## add a species column first half is GWfG as that come first in the alphabet
fitGWfG$Species <- "White-fronted"



#### Plot the output of the separate GBG and GWfG models on the same plot
fit2 <- rbind(fitGBG, fitGWfG)
## Now plot using ggplot
ggplot() + 
  geom_line(data=fit2, aes(x= dist_to_shot, y = fit, group = Species, colour = Species), size = 1.25)  +
  geom_ribbon(data = fit2, aes(x=dist_to_shot, ymin = lower, ymax = upper, group = Species, colour = Species), 
              alpha = 0.3, colour = NA, fill = "grey") + 
  ylab("Speed difference (m/s) [(t+1)-(t-1)]") + xlab("Distance to shot (m)") +
  theme_bw() +
  scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  xlim (0, 3000) +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=14,), 
        panel.grid.minor.x = element_blank())

setwd("~/Shooting disturbance/Plots for thesis chapter")
ggsave("Displacement as function of distance to shot [(t+1)-(t-1)].png", 
       width = 25, height = 15, units = "cm")








##                                                                           ##
#### 12. Use bootstrapped samples to estimate uncertainty in threshold points ####
##                                                                           ##

## If a I create bootstrapped samples and then recalculate the threshold point I should get
## an estimate for the uncertainty in the threshold point

## First set up values that control iterations and sampling for the bootstrap
n_sim <- 1000 # number of bootstrap samples to create and calculate the thresholds
dist_div <- 10 # increment in number of meters to search through to find threshold

## create a data frame that stores the best threshold for each sample

## create a data frame that stores the AICc values for each trialed threshold in each bootstrapped samples


## 
for(j in 1:n_sim){
  
  ## set seed so i get the same results if I run again
  set.seed(1212)
  
  ## draw a bootstrap sample from main data set
  
  
}




## might actually be easier to write a function that take my data, then finds the threshold and return the distance where the threshold point was found
## https://www.datacamp.com/community/tutorials/bootstrap-r
## https://data-flair.training/blogs/bootstrapping-in-r/
## https://towardsdatascience.com/a-practical-guide-to-bootstrap-with-r-examples-bd975ec6dcea

## boot package will be far more flexible and will probably good to learn how to use it for the future
## In machine learning, it is common to use a sample size that is the same as the original dataset.
## If the dataset is enormous and computational efficiency is an issue, smaller samples can be used, such as 50% or 80% of the size of the dataset.





## Function to calculate threshold point for a disturbance model
Threshold_finder <- function(data = data){
  
  ## create data frame to put the output in
  output_boot <- data.frame(model=as.numeric(),AIC=numeric())
  
  ## create sequence of distances, each of which will be used as the break point in the model
  dist_seq_boot <- seq(0, 3000, 10)
  
  ## loop through each threshold value individual and save model AICc in ouput data set
  for(i in 1:length(dist_seq_boot)) {
    
    ## run the model
    model <- lmer(speed_dif ~ threshold1(min_shot_dist, dist_seq_boot[i]) + TOD_sc + year + from_nov1_sc + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                  data= data,
                  REML=FALSE)
    
    output_boot[i,]<-c(dist_seq_boot[i], MuMIn::AICc(model))
    
  }
  
  ## Which model has the lowest AICc
  threshold_row <- as.numeric(which.min(output_boot$AIC))
  
  ## return the distance for the model with the lowest AICc
  Threshold_point <- output_boot$model[threshold_row]
  Threshold_point

}


## using the boot function calcualte a bootstrap CI for the threshold point
## need to think about on what levels i do the sampling
## I.e. just sample from all data, stratify so just sample within individuals, sampling unit as each individual
boot_guess <- boot::boot(data = GBG_pairs, statistic = Threshold_finder, R = 100, sim = "parametric")

summary(boot_guess)



##                                                                                                 ##
#### 13. Look at how the shooting response changes depending on which fix is chosen as the control ####
##                                                                                                 ##

## try this section for just GBG
GBG_allpairs <- filter(final_pairs4, species == "GBG")

## Run model, data set was already prepared in an earlier section
## first make sure the control column is set as a factor
GBG_allpairs$control <- as.factor(GBG_allpairs$control)

## In these models want to see how using the fix before or after the shooting event effect relationship between speed difs and distance to shot
## use the same 200m threshold we found earlier and drop days since nov 1st as it never does anything in the modes
Comp_mod <- lmer(speed_dif ~ threshold1(min_shot_dist, 400)*control + TOD_sc + year + (1|Shoot_ID1), 
                data = GBG_allpairs, 
                REML = F)

## check model summary
summary(Comp_mod) ## singlualr fit if put Tag ID in as a random effect...
drop1(Comp_mod, test = "Chi")

## Check model performance and assumptions
model_performance(Comp_mod) # marginal R squared is pretty poor
check_model(Comp_mod) # not that bad, some of the check indicate that the response is overdispersed

## Model Selection

## Run model selection
## change default "na.omit" to prevent models being fitted to different datasets
options(na.action = "na.fail") 

## create all candidate models using dredge, specify any dependencies (trace shows progress bar)
dredge_set <- MuMIn::dredge(Comp_mod, trace = 2)
nested_set <- subset(dredge_set, !MuMIn::nested(.), recalc.weights=T)
delta6_set <- subset(nested_set, delta<=6, recalc.weights=T)


## use the effects package, and ggplot to plot the model predictions
## first the effects for each predicitor
divisions <- 200
top_comp_effects <- predictorEffects(Comp_mod, focal.levels = divisions)
plot(top_comp_effects[1])

## now extract the fits for the first variable and bind them together
effectscomp <- top_comp_effects[1]
fitcomp <- as.data.frame(cbind(effectscomp[["min_shot_dist"]][["fit"]], effectscomp[["min_shot_dist"]][["lower"]], 
                            effectscomp[["min_shot_dist"]][["upper"]], effectscomp[["min_shot_dist"]][["x"]][["min_shot_dist"]]))
## change the names to something meaningful
setnames(fitcomp, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "dist_to_shot"))
## add a species column first half is GBG as that come first in the alphabet
fitcomp$Control <- NA
fitcomp$Control[1:divisions] <- "After disturbance"
fitcomp$Control[(divisions+1):(divisions*2)] <- "Before disturbance"

## Now plot using ggplot
ggplot() + 
  geom_line(data=fitcomp, aes(x= dist_to_shot, y = fit, group = Control, colour = Control), size = 1.25)  +
  geom_ribbon(data = fitcomp, aes(x=dist_to_shot, ymin = lower, ymax = upper, group = Control, colour = Control), 
              alpha = 0.275, colour = NA, fill = "grey") + 
  ylab("Speed difference (m/s)") + xlab("Distance to shot (m)") +
  theme_bw() +
  scale_colour_manual(values=c("#746AB0", "#FFCE30")) +
  xlim (0, 2000) +
  ggtitle("GBG response only") +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank())










































## OLD CODE GRAVEYARD ##



final_pairsGWfG <- filter(final_pairs4, species == "GBG")
## loop through each threshold value individual and save model AICc in ouput data set
for(i in 1:length(dist_seq)) {
  
  message(i, " out of ", length(dist_seq))
  
  ## In these models want to see how using the fix before or after the shooting event effect relationship between speed difs and distance to shot
  Comp_mod2 <- lmer(speed_dif ~ threshold1(min_shot_dist, dist_seq[i])*control + TOD_sc + year + from_nov1_sc + I(min_shot_dist_sc^2) + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                    data = final_pairsGWfG, 
                    REML = F)
  
  
  outputGWfG[i,]<-c(dist_seq[i], MuMIn::AICc(Comp_mod2))
  
}

outputGWfG
plot(outputGWfG$model, outputGWfG$AIC)
which.min(outputGWfG$AIC) ## row 10 which is 450m


#### Plot output of the top threshold model
top.model <- lmer(speed_dif ~ threshold1(min_shot_dist, 200)*control + TOD_sc + year + from_nov1_sc + I(min_shot_dist_sc^2) + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                  data = final_pairsGWfG, 
                  REML = F)

## use the effects package, easier than ggplot and saves times
top_mod_effects <- predictorEffects(top.model)
plot(top_mod_effects[1])


Pair_fix$Pair_ID <- 1:nrow(pairs2)
Dist_fix$Disturb <- 1
Pair_fix$Disturb <- 0

## now bind the two data sets on top of one another
## change one of the column names so the bind works
setnames(Pair_fix, old = "ID_prev", new = "ID")
final_pairs <- rbind(Dist_fix, Pair_fix)

## Finally join the full data for each of these fixes
final_pairs <- left_join(final_pairs, winter_gps_prox)



Ind_mod_top <- lmer(speed_dif ~ poly(Shoot_ID1_dist_to,2) + TOD_ + (1|Tag_ID) + (1|Shoot_ID1), 
                    data = final_pairs2, 
                    REML = F)

zz <- predictorEffects(Ind_mod_top)
plot(zz)
dd <- zz[["Shoot_ID1_dist_to"]][["data"]]
#create sequence of hour values
DistValues <- seq(0, 2000, 10)

#create list of predicted happines levels using quadratic model
bb <- predict(Ind_mod_top, aa)

aa <- list(Shoot_ID1_dist_to=DistValues, dist_shot_sq=DistValues^2)

#create scatterplot of original data values
plot(data$hours, data$happiness, pch=16)
#add predicted lines based on quadratic regression model
lines(hourValues, happinessPredict, col='blue')





#### OLD CODE TO RUN NON-THRESHOLD MODEL

## Run model
unique(final_pairs4$control)
Before_pairs <- filter(final_pairs4, control == "Shoot_T-1")
## In my mind this is modeling the speed differences between disturbed fixes and their controls
## as a function of distance from shot and other parameters....
Ind_mod <- lmer(speed_dif ~ min_shot_dist_sc*species + TOD_sc + year + from_nov1_sc + I(min_shot_dist_sc^2) + I(from_nov1_sc^2) + (1|Shoot_ID1), 
                data = Before_pairs, 
                REML = F)

## check model summary
summary(Ind_mod) ## singlualr fit if put Tag ID in as a random effect...
drop1(Ind_mod, test = "Chi")

## Check model performance and assumptions
model_performance(Ind_mod) # marginal R squared is pretty poor
check_model(Ind_mod) # not that bad, some of the check indicate that the response is overdispersed

## Model Selection

## Run model selection
## change default "na.omit" to prevent models being fitted to different datasets
options(na.action = "na.fail") 

## create all candidate models using dredge, specify any dependencies (trace shows progress bar)
dredge_set <- MuMIn::dredge(Ind_mod, trace = 2)
nested_set <- subset(dredge_set, !MuMIn::nested(.), recalc.weights=T)
delta6_set <- subset(nested_set, delta<=6, recalc.weights=T)



####


## Plot basically using just the effects package and base plotting
## Note I could probably extract the values in top_mod_effects and plot them in ggplot
Ind_mod_top <- lmer(speed_dif ~ poly(min_shot_dist, 2) + TOD_ + (1|Shoot_ID1), 
                    data = Before_pairs, 
                    REML = F)

top_mod_effects <- predictorEffects(Ind_mod_top)
plot(top_mod_effects)



## ***NOTE: this bit of code below does not plot the right relationship ***####
## Run the top model from model selection
## run with unsealed values so plot axis make sense
final_pairs2$dist_shot_sq <- (Before_pairs$Shoot_ID1_dist_to)^2
Ind_mod_top <- lmer(speed_dif ~ poly(Shoot_ID1_dist_to, 2) + dist_shot_sq + (1|Tag_ID) + (1|Shoot_ID1), 
                    data = Before_pairs, 
                    REML = F)


## extract the fitted values plus CIs using the effects package
effectz <- effects::effect(term= c("poly(Shoot_ID1_dist_to, 2)"), mod= Ind_mod_top, xlevels= 100)
effectz2 <- as.data.frame(effectz)

ggplot() + 
  geom_point(data=Before_pairs, aes(x=Shoot_ID1_dist_to, y=speed_dif, colour = "blue"), size =1) + 
  geom_line(data=effectz2, aes(x= dist_shot_sq, y = fit), size = 1.25)  +
  geom_ribbon(data = effectz2, aes(x=dist_shot_sq, ymin = lower, ymax = upper), alpha = 0.5, colour = NA, fill = "grey") + 
  ylab("Speed difference/m/s") + xlab("Distance to shot/m") +
  theme_bw() +
  ylim(3,-2)


