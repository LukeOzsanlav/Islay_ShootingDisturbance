## Luke Ozsanlav-Harris

## 13/01/2022

## Create a data set with all of the winter GPS fixes for all Scottish GWfG
## Where the birds were within a certain spatio temporal threshold to a shooting event is labeled in the data set
## Data set is for an undergraduate prject

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
library(raster)



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
setwd("~/Acc/ALL_ORNI_DATA/Cleaned CSVs")

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

##  order by TagID and by timestamp then add a unique ID to each fix
winter_gps <- winter_gps[order(winter_gps$Tag_ID, winter_gps$UTC_datetime),]
winter_gps$ID <- 1:nrow(winter_gps)





##                                                       ##
#### 1.4 Read out GWfG tag codes for associations matrix ####
##                                                       ##

## Read in tag info so I can add sexes to the tags
setwd("~/Additional data files")
tag_info <- fread("Tagged bird summary data.csv")

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
write_csv(codes2, file = "~/Shooting disturbance/List of tags used/GWfG_winter_codes.csv")







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
Fields <- st_read(dsn = "~/Shooting disturbance/88090_ISLAY_GMS_FIELD_BOUNDARY",
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
system.time(Prox_list <- Shoot_proximity(TD = winter_resamp, SD = winter_shots_cent2, time_tresh = 120, dist_thresh = 600, field_crs = Field_centres))

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





##                                    ##
#### 8. Add field ID to the GPS fixes ####
##                                    ##


## Convert the GPS data into a sf object
winter_gps_sf <- st_as_sf(winter_gps_prox, coords = c("Longitude", "Latitude"), 
                          crs = 4326, agr = "constant")

## now convert to the same crs as the field data
winter_gps_sf <- st_transform(winter_gps_sf, crs = st_crs(Fields))

## extract the polygon values to the points
Field_labels <- st_intersection(winter_gps_sf, Fields)

## Just extract the Field IDs
Field_labels <- subset(Field_labels, select = c("ID", "FIELD_ID"))
st_geometry(Field_labels) <- NULL

## Now join Field IDs back to the mina data set
winter_gps_fields <- full_join(winter_gps_sf, Field_labels, by = "ID")


## Plotting function if needed
ggplot() +
  geom_sf(data = Fields, fill="lightgrey") + # coast3 is just shapefile of the worlds coastline, use Islay field boundaries instead?
  geom_sf(data = winter_gps_100, fill="black") + # use some other data set here, maybe field centroids or shooting locations
  theme_light() + # set a theme
  theme(panel.grid.major = element_blank()) # remove the guidelines


## read out the data set
#write_csv(winter_gps_fields, file = "~/Shooting disturbance/Datasets for Marta/WinterGPS_forMarta.csv")





##                                        ##
#### 9. Add habitat data to the GPS fixes ####
##                                        ##


##                                    ##
#### 9.1 Read in the WWT habitat data ####
##                                    ##

## This data is from May 2013
## Read in the csv with the WWT habitat data
setwd("~/Shooting disturbance/Habitat data")
WWTcsv <- fread("WWT habitat data 2014.csv")

## select just the column i want from the WWT habitat data and rename Field ID
WWTHabitat <- WWTcsv[ ,c("field_id", "domveg")]
setnames(WWTHabitat, old = "field_id", new = "FIELD_ID")

##Check for duplicated observations in WWT habitat data
Dups <- WWTHabitat %>% 
        dplyr::select(FIELD_ID) %>%
        duplicated()
sum(Dups) 
WWTHabitat <- WWTHabitat %>% filter(!Dups)


## bind together the Field sf object with the vegetation categories from WWT
Fields_WWTHab <- left_join(Fields, WWTHabitat, by = "FIELD_ID")
summary(duplicated(Fields_WWTHab$FIELD_ID)) # check for duplicates again

## plot this habitat data
ggplot() +
  geom_sf(data = Fields_WWTHab, aes(fill=domveg)) + 
  theme_light() + # set a theme
  theme(panel.grid.major = element_blank())

## Make a list of the WWT habitat classes
WWTClasses <- unique(WWTHabitat$domveg)






##                                      ##
#### 9.2 Read in the Aimee habitat data ####
##                                      ##


## This data is from 2019
## read in the raster with the habitat data
Scot_hab <- raster("~/Shooting disturbance/Habitat data/ukregion-scotland.tif")

## Crop extent to just Islay
islay <- extent(114500, 150000, 639000,680000) # This extent applies for CEH land use data and EUNIS
Islay_hab <- crop(Scot_hab, islay)
plot(Islay_hab)

## create a list of unique values
IslayClasses <- as.data.frame(unique(values(Islay_hab)))

Islay_hab2 <- Islay_hab
Islay_hab2[values(Islay_hab)<7] <- 0
Islay_hab2[values(Islay_hab)>7] <- 0
plot(Islay_hab2)





##                                              ##
#### 9.3 Convert WWT habitat data into a Raster ####
##                                              ##

## create empty raster with same dimensions as Islay_hab data
empty_rast <- Islay_hab
values(empty_rast) <- NA

## Make sure the WWT habitat polygons have the same crs as the empty raster
Fields_WWTHabTr <- st_transform(Fields_WWTHab, crs = st_crs(empty_rast))

## Give each habitat a numeric value
domveg <- unique(WWTHabitat$domveg)
HabCodes <- as.data.frame(domveg)
HabCodes$HabCode <- 1:nrow(HabCodes)
Fields_WWTHabTr <- full_join(Fields_WWTHabTr, HabCodes, by = "domveg")
Fields_WWTHabTr$HabCode <- as.numeric(Fields_WWTHabTr$HabCode)

## No convert the WWT habitat data to a raster
## NOTE: Need to look into the FUN argument in this function, may chnage results slightly where a pixel is overlapped by two polygons
WWT_HabRaster <- rasterize(Fields_WWTHabTr, y = empty_rast, field = "HabCode")
plot(WWT_HabRaster)





##                           ##
#### 10. Split tracking data ####
##                           ##

## Split the trucking data into two halves depending on whether they are closer to May 2013 or June 2019
## The split will be in Summer 2016 

Pre2016GPS <- filter(winter_gps_fields, UTC_datetime < "2016-06-15 23:00:00")
Post2016GPS <- filter(winter_gps_fields, UTC_datetime > "2016-06-15 23:00:00")





##                                                   ##
#### 11. Assign habitats to pre Summer 2016 GPS data ####
##                                                   ##

## Make sure the GPS fixes have the same crs as the Islay hab raster (whihc is the sama as the WWT raster)
Pre2016GPS <- st_transform(Pre2016GPS, crs = st_crs(Islay_hab))

## extract the raster values to the points for both habitat data sets
Pre2016GPS$HabitatCodeWWT <- raster::extract(x = WWT_HabRaster, y = Pre2016GPS, method = "simple")
Pre2016GPS$HabitatCodeIsla <- raster::extract(x = Islay_hab, y = Pre2016GPS, method = "simple")





##                                                    ##
#### 12. Assign habitats to post Summer 2016 GPS data ####
##                                                    ##

## Make sure the GPS fixes have the same crs as the Islay hab raster
Post2016GPS <- st_transform(Post2016GPS, crs = st_crs(Islay_hab))

## extract the raster values to the points
Post2016GPS$HabitatCode <- raster::extract(x = Islay_hab, y = Post2016GPS, method = "simple")





##                                                        ##
#### 13. Convert habitat codes and join two GPS data sets ####
##                                                        ##


## Read in the two conversion tables
setwd("~/Shooting disturbance/Habitat data")
WWTConv <- read.csv("HabitatConversion_WWTdata.csv", header = T); WWTConv <- WWTConv[,1:2]
IslaConv <- fread("HabitatConversion_Aimeedata.csv")

## Need to add the numeric codes to the WWT conversion
WWTConv <- full_join(WWTConv, HabCodes, by = "domveg")



## Add the converted habitat codes to the post 2016 GPS habitat data
IslaConv <- subset(IslaConv, select = c("RasterValue", "HabClass"))
setnames(IslaConv, old = "RasterValue", new = "HabitatCode")
Post2016GPS2 <- left_join(Post2016GPS, IslaConv, by = "HabitatCode")



## Add the habitat codes to the WWT habitat data
## Have to add the codes for the WWT data dn the Islay data
## First add the converted habitat codes for the HabitatCodeIsla column
colnames(Pre2016GPS)
setnames(IslaConv, old = "HabClass", new = "HabClassIsla")
Pre2016GPS2 <- left_join(Pre2016GPS, IslaConv,  by = c("HabitatCodeIsla" = "HabitatCode"))

## Now add the converted habitat codes for the HabitatCodeWWT column
colnames(WWTConv)
WWTConv <- subset(WWTConv, select = c("HabClass", "HabCode"))
setnames(WWTConv, old = "HabClass", new = "HabClassWWT")
Pre2016GPS2 <- left_join(Pre2016GPS2, WWTConv,  by = c("HabitatCodeWWT" = "HabCode"))


## Finally if there isn't a habitat for the WWT data then just use the one in the Islay data set
#Pre2016GPS2$HabClass <- ifelse(is.na(Pre2016GPS2$HabClassWWT)==T & is.na(Pre2016GPS2$HabClassIsla)==F, Pre2016GPS2$HabClassIsla, Pre2016GPS2$HabClassWWT)
Pre2016GPS2$HabClass <- Pre2016GPS2$HabClassIsla



## There are mismatches  between the two habitat data sets as one is at the pixel level and the other is at the field level
## For Martas stuff I am just going to send the Islay habitat data from 2019


## Bind together the two data sets
colnames(Pre2016GPS2); colnames(Post2016GPS2)
Pre2016GPS2$HabitatCodeWWT <- NULL
Pre2016GPS2$HabitatCodeIsla <- NULL
Pre2016GPS2$HabClassIsla <- NULL
Pre2016GPS2$HabClassWWT <- NULL
Post2016GPS2$HabitatCode <- NULL

AllGPS_withHabitat <- rbind(Pre2016GPS2, Post2016GPS2)
write_csv(AllGPS_withHabitat, file = "~/Shooting disturbance/Datasets for Marta/WinterGPS_forMarta.csv")
unique(AllGPS_withHabitat$HabClass)
