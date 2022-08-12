## Luke Ozsanlav-Harris

## Identify disturbance events in tracking data that are consequence of shooting
## Find fixes that are spatio temporally associated with shooting events
## Determine whether the birds were displaced during the gap when the shooting event occurred

## Fixes:
## Try initially picking the fixes nearest the shooting event
## 



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
library(svMisc)



##                                          ##
#### 1.  Read in and clean Ecotone GPS data ####
##                                          ##


## set working directory
setwd("~/Ecotone Data/Cleaned Data file")

## read in data
## Will need to read in allthe data at a later point but just start with one file for now
#orn_data <-list.files(pattern = "*.csv") %>% map_df(~fread(.))
eco_data <- fread("Ecotone_all_clean.csv")


#### Data Cleaning ####
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

#### Data Cleaning ####
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

##  Add row ID to the tracking data
winter_gps$ID <- 1:nrow(winter_gps)







##                                                 ##
#### 2. Read in shooting data and field boundaries ####
##                                                 ##

## read in the shooting logs
setwd("~/Shooting disturbance")
Logs <- fread("All_logs_cleaned.csv")

## create a timestamp column for the shooting logs 
## think the ones that fail to parse are the ones were the time is missing and it is just a date
Logs$timestamp <- ymd_hms(paste0(Logs$Date_Cl, " ", Logs$Time_Cl))


## Now filter the logs for the same time period as the shooting data
winter_logs <- filter(Logs, timestamp > min(winter_gps$UTC_datetime) & timestamp < max(winter_gps$UTC_datetime))

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






##                                    ##
#### 3. Manipulate spatial field data ####
##                                    ##


## calculate the centroid of each field
Field_centres <- st_centroid(Fields)

## now create a data set with the shooting events but where each one has the buffer for that field
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
#### 4. Label fixes associated with each shooting event ####
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
    for (i in 1:nrid) {
      
      progress(i, max.value = nrid)
      
      
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
system.time(Prox_list <- Shoot_proximity(TD = winter_gps, SD = winter_shots_cent2, time_tresh = 120, dist_thresh = 600, field_crs = Field_centres))

## bind the lists together into one data frame
All_prox <- plyr::ldply(Prox_list)







##                                                            ##
#### 5. Deal with multi associations and add shooting times ####
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


##** UP TO HERE **##
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
#### 6. Add shooting associations back to the full tracking data set ####
##                                                                   ##


## Join the proximity data set back to the full gps set
colnames(Shoot_prox)
Shoot_prox2 <- subset(as.data.frame(Shoot_prox), select = -c(Tag_ID, UTC_datetime))
winter_gps_prox <- full_join(winter_gps, Shoot_prox2, by = "ID")

## assign all NAs in the spatial overlap column a 0
winter_gps_prox$spatial_overlap <- ifelse(is.na(winter_gps_prox$spatial_overlap) == T, 0, winter_gps_prox$spatial_overlap)


#### Calculate step lengths between the fixes ####

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

## Now determine which fixes should be in an exposure set
min_time <- 120 # should be in mins
winter_gps_prox$Expo_set <- ifelse(winter_gps_prox$spatial_overlap == 1 & winter_gps_prox$time_dif <= min_time, 
                                   1, 0)

## write out the data set where the GPS data is labeled with the shooting events
#winter_gps_prox$Expo_set <- NULL
#write_csv(winter_gps_prox, file = "~/Shooting disturbance/WinterGPS_labeled_with_shooting.csv")







##                                                             ##
#### 7. Identify likely displacement events in tracking data ####
##                                                             ##


## Use a combination of step length and change in field code to identify displacement events
## labeling GPS with field codes may be slow so may want to filter data to bird days that were in proximity to a shooting event


#### Subset tracking data to tag days with shooting events ####

## extract just the bird days were there was proximity to a shooting event (These are in the Shoot Prox data set)
## add bird day column to the full GPS data set and the shooting proximity data sets
Shoot_prox$tag_day <- paste0(Shoot_prox$Tag_ID, "_", as.Date(Shoot_prox$UTC_datetime))
winter_gps_prox$tag_day <- paste0(winter_gps_prox$Tag_ID, "_", as.Date(winter_gps_prox$UTC_datetime))

## create a list of the unique tag days
unique_days <- unique(Shoot_prox$tag_day)

## Now filter out the days were birds were in proximity to shooting events
gps_shoot_days <- filter(winter_gps_prox, tag_day %in% c(unique_days))




#### Label fixes with field codes ####
  
## convert these GPS fixes to an sf object
gps_shoot_days_sf = st_as_sf(gps_shoot_days, coords = c("Longitude", "Latitude"), 
                         crs = 4326, agr = "constant")

## now convert to the same crs as the field data
gps_shoot_days_sf <- st_transform(gps_shoot_days_sf, crs = st_crs(Fields))

## extract which fields fixes lie in
## if they don't lie in a field then they are left out fo the final data set
overlap_ <- st_intersection(gps_shoot_days_sf, Fields)

## select just the GPS ID column and the field code
overlap_sub <- overlap_ %>% 
                  as.data.frame() %>% 
                  subset(select = c("ID", "FIELD_ID"))

## now join this back onto the full gps data set
gps_shoot_days2 <- left_join(gps_shoot_days, overlap_sub, by = "ID")





##                                                                                  ##
#### 8. Find displacement events with a shooting event in the inter-fix intervals ####
##                                                                                  ##

## Want to find the fix immediately before the shooting event

## First find the particular inter fix period were the shooting event took place
## Then see if the bird traveled more than x metes and changed field code

## Set thresholds for displacement
Dist_Thresh <- 0.5 # threshold in meters of a movement to be a displacement
Time_thresh <- 120 # threshold in minutes for a inter fix interval to be too large to discount displacement


## Split up the data by which slot the shooting data is in
Slot1 <- gps_shoot_days2[is.na(gps_shoot_days2$Shoot_ID1) == F,]
Slot2 <- gps_shoot_days2[is.na(gps_shoot_days2$Shoot_ID2) == F,]
Slot3 <- gps_shoot_days2[is.na(gps_shoot_days2$Shoot_ID3) == F,]
Slot4 <- gps_shoot_days2[is.na(gps_shoot_days2$Shoot_ID4) == F,]


#### Label events where birds have likely been displaced ##
#### Shooting slot 1 
## apply distance and time rules to label
colnames(Slot1)
Slot1$Shoot_displac <- ifelse(Slot1$Step_len > Dist_Thresh &
                              Slot1$time_dif < Time_thresh &
                              is.na(Slot1$Shoot_ID1) == F, 
                              1, 0)

#Slot1_1 <- subset(Slot1,select = c("tag_year", "Tag_ID", "Step_len", "time_dif", "tag_day", "FIELD_ID", "Shoot_displac", "UTC_datetime", "Latitude", "Longitude", "year_day", "ID") )


#### Shooting slot 2 
## apply distance and time rules to label
Slot2$Shoot_displac <- ifelse(Slot2$Step_len > Dist_Thresh &
                              Slot2$time_dif < Time_thresh &
                              is.na(Slot2$Shoot_ID2) == F, 
                              1, 0)

#### Shooting slot 3 
## apply distance and time rules to label
Slot3$Shoot_displac <- ifelse(Slot3$Step_len > Dist_Thresh &
                              Slot3$time_dif < Time_thresh &
                              is.na(Slot3$Shoot_ID3) == F, 
                              1, 0)

#### Shooting slot 4 
## apply distance and time rules to label
Slot4$Shoot_displac <- ifelse(Slot4$Step_len > Dist_Thresh &
                              Slot4$time_dif < Time_thresh &
                              is.na(Slot4$Shoot_ID4) == F, 
                              1, 0)


## bind the data back together form all of the slots
Shoot_exposure <- rbind(Slot1, Slot2, Slot3)
table(Slot2$Shoot_displac)




##                                                                     ##
#### 9. Model how distance from shot effect displacement probability ####
##                                                                     ##

## just going to do it with the first shooting slot but need to expand it out to all four
## I think ultimately it will be easier to just repair the fixes that are associated with multiple shooing events in the data set instead of multiple columns

## Run a GLMM to model probability of disturbance


## filter out just the set considered for exposure, some were ruled out as the gap between fixes was too long
Slot1 <- filter(Slot1, Expo_set == 1)

## convert logit to probability
logit2prob <- function(logit){
                odds <- exp(logit)
                prob <- odds / (1 + odds)
                return(prob)
}

## Make sure variables are the right class 
Slot1$Shoot_ID1_dist_to <- as.numeric(Slot1$Shoot_ID1_dist_to)
Slot1$Shoot_ID1_shots <- as.numeric(Slot1$Shoot_ID1_shots)
Slot1$Shoot_ID1_time_dif <- as.numeric(Slot1$Shoot_ID1_time_dif)

## Run the GLMM
## add time of day, time of year, year
colnames(Slot1)
table(Slot1$Tag_ID)

## add additional columns to the model
## year column
Slot1$year <- year(Slot1$UTC_datetime)
## time of day column
Slot1$TOD_ <- as.numeric(hour(Slot1$Shoot_ID1_time)) + as.numeric(minute(Slot1$Shoot_ID1_time)/60)
## days Since November 1st column
Slot1$Since_Nov1 <- ifelse(Slot1$year_day > 250, 
                                   Slot1$year_day - 305, 
                                   Slot1$year_day + 60)




## scale continuous variables
Slot1$dist_to_shot_sc <- scale(Slot1$Shoot_ID1_dist_to)
Slot1$Since_Nov1_sc <- scale(Slot1$Since_Nov1)
Slot1$TOD_sc <- scale(Slot1$TOD_)
Slot1$Shots_fired_sc <- scale(Slot1$Shoot_ID1_shots)
Slot1$time_from_shot_sc <- scale(Slot1$Shoot_ID1_time_dif)
Slot1$time_dif_sc <- scale(Slot1$time_dif)

## remove an NA rows 
Slot1_1 <- Slot1 %>% drop_na(Shoot_displac)
table(Slot1_1$Shoot_displac) ## size of data set

Diplacement_mod <- glmer(Shoot_displac ~ dist_to_shot_sc + time_dif_sc + Since_Nov1_sc + TOD_sc + Shots_fired_sc +
                                          (1|Tag_ID),
                         family = binomial(link=logit),
                         data = Slot1_1, 
                         control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e8)))

summary(Diplacement_mod)


## check the model with DHARMa
simulationOutput <- simulateResiduals(fittedModel = Diplacement_mod, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

## plot residuals vs predictors
plotResiduals(simulationOutput, Slot1_1$dist_to_shot_sc) ## bir problematic but might just be alright
plotResiduals(simulationOutput, Slot1_1$Since_Nov1_sc)
plotResiduals(simulationOutput, Slot1_1$TOD_sc)
plotResiduals(simulationOutput, Slot1_1$time_dif_sc)
plotResiduals(simulationOutput, Slot1_1$Shots_fired_sc)


## Run model selection using MuMIN
## change default "na.omit" to prevent models being fitted to different datasets
options(na.action = "na.fail") 

## create all candidate models using dredge, specify any dependencies (trace shows progress bar)
ms2 <- MuMIn::dredge(Diplacement_mod, trace = 2)
ms2_sub <- subset(ms2, !MuMIn::nested(.), recalc.weights=T)

## run the top model
top_mod <- glmer(Shoot_displac ~ dist_to_shot_sc + time_dif_sc + Shots_fired_sc +
                                 (1|Tag_ID),
                family = binomial(link=logit),
                data = Slot1_1, 
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e8)))
summary(top_mod)
drop1(top_mod, test = "Chisq")
r.squaredGLMM(top_mod); confint(top_mod)

## Marginal R2 represents the variance explained by the fixed effects
## conditional R2 is the variance explained by the entire model (i.e. the fixed and random effects). 
## As a consequence, the marginal R2 cannot be higher than the conditional R2.


## plot the effects from the top model
## extract the fitted values plus CIs using the effects package
effectz <- effects::effect(term= "dist_to_shot_sc", mod= top_mod, xlevels= 100)
effectz2 <- as.data.frame(effectz)


## create the plot
ggplot() + 
  geom_point(data=Slot1_1, aes(x=dist_to_shot_sc, y=Shoot_displac), size =1.5) + 
  geom_line(data=effectz2, aes(x= dist_to_shot_sc, y = fit), size = 1.25)  +
  geom_ribbon(data = effectz2, aes(x=dist_to_shot_sc, ymin = lower, ymax = upper), alpha = 0.5, colour = NA, fill = "grey") + 
  ylab("Probability of displacement") + xlab("Distance from shooting/km") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        panel.grid.minor.x = element_blank(), strip.text.x = element_text(size = 12, face = "bold"))


## plot the effects from the top model
## extract the fitted values plus CIs using the effects package
effectz3 <- effects::effect(term= "time_dif_sc", mod= top_mod, xlevels= 100)
effectz4 <- as.data.frame(effectz3)

## create the plot
ggplot() + 
  geom_point(data=Slot1_1, aes(x=time_dif_sc, y=Shoot_displac), size =1.5) + 
  geom_line(data=effectz4, aes(x= time_dif_sc, y = fit), size = 1.25)  +
  geom_ribbon(data = effectz4, aes(x=time_dif_sc, ymin = lower, ymax = upper), alpha = 0.5, colour = NA, fill = "grey") + 
  ylab("Probability of displacement") + xlab("Inter fix interval") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        panel.grid.minor.x = element_blank(), strip.text.x = element_text(size = 12, face = "bold"))



## plot the effects from the top model
## extract the fitted values plus CIs using the effects package
effectz5 <- effects::effect(term= "Shots_fired_sc", mod= top_mod, xlevels= 100)
effectz6 <- as.data.frame(effectz5)

## create the plot
ggplot() + 
  geom_point(data=Slot1_1, aes(x=Shots_fired_sc, y=Shoot_displac), size =1.5) + 
  geom_line(data=effectz6, aes(x= Shots_fired_sc, y = fit), size = 1.25)  +
  geom_ribbon(data = effectz6, aes(x=Shots_fired_sc, ymin = lower, ymax = upper), alpha = 0.5, colour = NA, fill = "grey") + 
  ylab("Probability of displacement") + xlab("Number of shots fired") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.text=element_text(size=12, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        panel.grid.minor.x = element_blank(), strip.text.x = element_text(size = 12, face = "bold"))










##















Time_filter_function <- function(x){
  
  gps_sub <- filter(winter_gps, UTC_datetime > (x$timestamp[j]-hours(1)) & UTC_datetime < (x$timestamp[j]+hours(1)))
  
  return(gps_sub)
  
}


Shooting_prox_function <- function(x){
  
  
  ## extract fixes 1 hour either side of the shooting event
  gps_sub <- filter(winter_gps, UTC_datetime > (x$timestamp[j]-hours(1)) & UTC_datetime < (x$timestamp[j]+hours(1)))
  
  if(nrow(gps_sub)==0){print("no fixes in time buffer")}
  else{
    ## convert these GPS fixes to an sf object
    gps_sub_sf = st_as_sf(gps_sub, coords = c("Longitude", "Latitude"), 
                          crs = 4326, agr = "constant")
    
    ## now convert to the same crs as the field data
    gps_sub_sf <- st_transform(gps_sub_sf, crs = st_crs(x))
    
    ## See if any of the fixes overlap with the field buffer
    Fix_overlap <- st_within(gps_sub_sf, x)
    
    ## assign a 1 if the point were inside the field buffer
    gps_sub_sf$spatial_overlap <- ifelse(as.matrix(Fix_overlap) == FALSE, 0, 1)
    
    return(gps_sub_sf)
  }
  
}



## extract the first 10 rows of the shootng buffer data for a practice with
winter_shots_sub <- winter_shots_buf2[1:100,]

system.time(Proximity <- winter_shots_sub %>% 
              group_by(ID) %>% 
              nest() %>% 
              mutate(Birds_near =map(data, Shooting_prox_function)) )


system.time(Proximity <- winter_shots_sub %>% 
             group_by(ID) %>% 
             nest() %>% 
             mutate(Birds_near =map(data, Time_filter_function)) )




