## Luke Ozsanlav-Harris
## 21/04/22

## Compare behavior pre and post shooting using classified acc data from Ornitela tags

## Packages required
library(tidyverse)
library(lubridate)
library(data.table)
library(zoo)
library(svMisc)

library(geosphere)
library(sf)
library(amt)

library(lme4)
library(DHARMa)
library(MuMIn)
library(effects)
library(performance)
library(glmmTMB)

## Result: Didn't found a difference in any behavior between the hour before and after shooting exposure
##         I tried varying the distance from the shooting event that birds had to be to be included in the data

## Problems: models aren't working the best, not 100% how to model this


##
#### 1. Read in data sets ####
##


##
#### 1.1 Read in classified ACC data ####
##
## Data Cleaning, extracting just the GPS data
orn_behav <- readRDS("tracking data/Classified_ACC.RDS")
table(orn_behav$device_id)

## read in meta data so i can just keep the Scotland birds
Meta <- fread("MetaData/Tagged bird summary data new.csv")
Meta <- filter(Meta, Ringing.location %in% c("ISLA", "WEST", "LOKE", "DYFI", "KINT") | S.N == 17795)

## filter the gps data for the device ids left
orn_behav <-filter(orn_behav, device_id %in% unique(Meta$S.N))

## parse timestamp
orn_behav$UTC_datetime <- ymd_hms(orn_behav$UTC_datetime)

##Check for duplicated observations (ones with same lat, long, timestamp, and device ID)
ind <- orn_behav %>% 
       dplyr::select(UTC_datetime, device_id) %>%
       duplicated()
sum(ind) 
orn_behav <- orn_behav %>% filter(!ind)


## Streamline the dataset for a first pass, first just have subset of the columns
setnames(orn_behav, old = "device_id", new = "Tag_ID")

## Now filter for the winter, create a year day column first
orn_behav$year_day <- yday(orn_behav$UTC_datetime)
winter_orn_behav <- filter(orn_behav, year_day >= 305 | year_day <= 92)
table(winter_orn_behav$Tag_ID)




##
#### 1.2 Read in Islay winter GPS data ####
##

GPS_data <- readRDS("Derived data/All_winter_GPS_with_habitat.RDS")

## only use the GPS data from birds with ACC data
ACC_tags <- unique(winter_orn_behav$Tag_ID)
GPS_set <- filter(GPS_data, Tag_ID %in% ACC_tags)

## give each fix its own unique ID number
GPS_set$ID <- 1:nrow(GPS_set)





##                                                 
#### 1.3 Read in shooting data and field boundaries ####
##                                                 


## read in the shooting logs
Logs <- fread("Shooting logs/All_logs_cleaned.csv")

## create a timestamp column for the shooting logs 
## think the ones that fail to parse are the ones were the time is missing and it is just a date
Logs$timestamp <- ymd_hms(paste0(Logs$Date_Cl, " ", Logs$Time_Cl))

## Now filter the logs for the same time period as the tracking data
winter_logs <- filter(Logs, timestamp > min(GPS_data$UTC_datetime) & timestamp < max(GPS_data$UTC_datetime))

## Now filter out just the shooting events
winter_shots <- filter(winter_logs, Shots_fired_Cl > 0)



## read in the shapefiles of the Islay field boundaries
Fields <- st_read(dsn = "Landcover Data/88090_ISLAY_GMS_FIELD_BOUNDARY",
                  layer = "AS_ISLAY_GMS_FIELD_BOUNDARY")
st_crs(Fields)

## remove any duplicates in the field data
ind3 <- duplicated(Fields$FIELD_ID)
sum(ind3) 
Fields <- Fields %>% filter(!ind3)




##                                                   
#### 2. Combine shooting data and spatial field data ####
##                                                   


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





##                                                   
#### 3. Label fixes associated with each shooting event ####
##                                                     


## Define function that converts timestamp immediately before shooting event the smallest value
## y needs to be the timestamps from the tracking data and x is just a single value from the shooting data sequence
difftime2 <- function(x=x, y=y){
  
  time_difs <- difftime(y, x, units='mins')*-1
  time_difs[time_difs<0] <- 1000000000
  return(as.vector(time_difs))
}

## Define function that will loop through each tag and find the nearest fix before each shooting event while the tag was active
## The fix prior to the shooting event has to be within a certain spatial and temporal buffer to the shooting event
## TD = tracking data. Need tag ID (Tag_ID), lat/long (Longitude/Latitude), timestamp (UTC_datetime) and unique ID for each for row (ID)
## SD = the shooting data. Need the time of the shot (timestamp), the field center (as an sf geometry) and unique ID for each for row (shot_ID)
## time_thresh = the max time in minutes that a fix can be before a shooting event
## dist_thresh = the max distance in meters the fix before a shooting event can be from the location of the shooting event
## field_crs = an sf object that has the crs of the field centers in SD

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
    
    if(length(shootseq)==0){next}
    
    
    ## FINDING NEAREST PRIOR FIX ##
    # find the first fix before the shooting event
    idx = sapply(shootseq, function(x) which.min(difftime2(y = Dtemp$UTC_datetime, x = x)))
    
    ## combine these fixes immediately before the shooting event with the shooting data
    aa <- Dtemp[idx,]
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
system.time(Prox_list <- Shoot_proximity(TD = GPS_set, SD = winter_shots_cent2, time_tresh = 60, dist_thresh = 2000, field_crs = Field_centres))

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

## Shoot prox now contains all of the GPS fixes that were in spatiotemporal proximity to shooting events
## This includes the time of the shooting event and how close the bird was to the shooting event




##
#### 7. Extract ACC behaviors before and after shooting events ####
##

## loop might work best
## extract all burst 1 hour before and after the shooting event then sum up the different behaviors
## Need unique pair IDs and the distance from shot
## later will need to add sex 


## Start by just using the instacnes were burds were associated with a single shooting event
Shoot_prox1 <- filter(Shoot_prox, is.na(Shoot_ID2_time)==T)

## create list of unique tag IDs
Tags <- unique(Shoot_prox1$Tag_ID)


for(i in 1:length(Tags)){
  
  message(i, " out of ", length(Tags)) 
  
  TagX <- filter(orn_behav, Tag_ID == Tags[i])
  Shoot_proxX <- filter(Shoot_prox1, Tag_ID == Tags[i])
  
  for(j in 1:nrow(Shoot_proxX)){
    
    progress(j, max.value = nrow(Shoot_proxX)) 
    
    ## filter out the data for the hour before the shooting event
    Sub_pre <- filter(TagX, UTC_datetime < Shoot_proxX$Shoot_ID1_time[j] & UTC_datetime > (Shoot_proxX$Shoot_ID1_time[j]-3600))
    Sub_pre$pair_ID <- j
    Sub_pre$pre_post <- "pre"
    Sub_pre$shoot_ID <- Shoot_proxX$Shoot_ID1[j]
    
    ## filter out the data for the hour before the shooting event
    Sub_post <- filter(TagX, UTC_datetime > Shoot_proxX$Shoot_ID1_time[j] & UTC_datetime < (Shoot_proxX$Shoot_ID1_time[j]+3600))
    Sub_post$pair_ID <- j
    Sub_post$pre_post <- "post"
    Sub_post$shoot_ID <- Shoot_proxX$Shoot_ID1[j]
    
    ## bind pre and post together
    sub <- rbind(Sub_pre, Sub_post)
    
    ## bind sub sections back together
    if(j==1){All_subs <- sub}else{All_subs <- rbind(All_subs, sub)}
    
  }
  
  if(i==1){All_tags <- All_subs}else{All_tags <- rbind(All_tags, All_subs)}
  
}





##
#### 8. Summarize behaviors before and after shooting event ####
##

## Create binary columns for each behaviour
All_tags$Flight <- ifelse(All_tags$Behaviour == "flight", 1, 0)
All_tags$Walk <- ifelse(All_tags$Behaviour == "walk", 1, 0)
All_tags$Graze <- ifelse(All_tags$Behaviour == "graze", 1, 0)
All_tags$Stat <- ifelse(All_tags$Behaviour == "stationary", 1, 0)
All_tags$Active <- ifelse(All_tags$Behaviour == "graze" | All_tags$Behaviour == "walk", 1, 0)

## Create summary of behaviors before and after shooting
colnames(All_tags)
BehSum <- All_tags %>% 
          group_by(Tag_ID, pair_ID, pre_post, shoot_ID) %>% 
          summarise(Flight = sum(Flight), 
                    Walk = sum(Walk),
                    Graze = sum(Graze),
                    Stat = sum(Stat),
                    Active = sum(Active),
                    Total = length(Behaviour))

## Add on the time of the shot and the distance from the shot
Shoot_sub <- subset(Shoot_prox1, select = c("Tag_ID", "Shoot_ID1", "Shoot_ID1_time", "Shoot_ID1_shots", "Shoot_ID1_dist_to"))
BehSum$Tag_ID <- as.character(BehSum$Tag_ID)
BehSum2 <- left_join(BehSum, Shoot_sub, by = c("shoot_ID"="Shoot_ID1", "Tag_ID"))




##
#### 9. Run model to compare behavior before and after shooting ####
##

## compare behavior 1 hour before and 1 hour after experiencing a shooting event

## filter out paired samples were there was less then 8 observations in the hour before or after shooting
BehSum3 <- filter(BehSum2, Total >=8)
BehSum3$Tag_pair <- paste0(BehSum3$Tag_ID, "_", BehSum3$pair_ID)
tag_pair_table <- as.data.frame(table(BehSum3$Tag_pair))
tag_pair_table <- filter(tag_pair_table, Freq == 2)
BehSum4 <- filter(BehSum3, Tag_pair %in% (tag_pair_table$Var1))


## Filter out pairs that were within a certain distance of the shootng event (tried 500m, 1000m and 2000m)
BehSum5 <- filter(BehSum4, Shoot_ID1_dist_to < 1000)


## initially tried a binoamial model (below) but had lots of problems with convergence
# Flight_mod <- glmer(cbind(Flight, Total) ~ pre_post+ (1|Tag_ID:pair_ID),
#                     data = BehSum5,
#                     family = binomial)






## Then tried to model it as a count, I.e. then number of graze instances in the hour sample
## I used an offset to account for the fact that the number of ACC samples differed between 8-10
## See this link for explanation of offset: https://stats.stackexchange.com/questions/11182/when-to-use-an-offset-in-a-poisson-regression

## I modelled counts of each of the 4 behaviours
## I havent controlled for things like Sex or day of year because of the paired nature of the smapling procedure. 

## Stationary model
Stat_mod <- glmmTMB(Stat ~ pre_post + offset(Total) + (1|Tag_ID:pair_ID),
                    data = BehSum5,
                    family = nbinom2,
                    control = glmmTMBControl(optimizer=optim,
                                             optArgs=list(method="BFGS", iter.max=1e5,eval.max=1e5)))
summary(Stat_mod)
drop1(Stat_mod, test = "Chi")
model_performance(Stat_mod) 
check_model(Stat_mod) 
check_overdispersion(Stat_mod)


## Flight model
Flight_mod <- glmmTMB(Flight ~ pre_post + offset(Total) + (1|Tag_ID:pair_ID),
                      data = BehSum5,
                      family = nbinom2,
                    control = glmmTMBControl(optimizer=optim,
                                             optArgs=list(method="BFGS", iter.max=1e5,eval.max=1e5)))
summary(Flight_mod)
drop1(Flight_mod, test = "Chi")
model_performance(Flight_mod) 
check_model(Flight_mod) 



##Walking model
Walk_mod <- glmmTMB(Walk ~ pre_post + offset(Total) + (1|Tag_ID:pair_ID),
                    data = BehSum5,
                    family = nbinom2,
                    control = glmmTMBControl(optimizer=optim,
                                               optArgs=list(method="BFGS", iter.max=1e5,eval.max=1e5)))
summary(Walk_mod)
drop1(Walk_mod, test = "Chi")
model_performance(Walk_mod) 
check_model(Walk_mod) 


## Grazing model
Graze_mod <- glmmTMB(Graze ~ pre_post + offset(Total) + (1|Tag_ID:pair_ID),
                    data = BehSum5,
                    family = nbinom2,
                    control = glmmTMBControl(optimizer=optim,
                                             optArgs=list(method="BFGS", iter.max=1e5,eval.max=1e5)))
summary(Graze_mod)
drop1(Graze_mod, test = "Chi")
model_performance(Graze_mod) 
check_model(Graze_mod) 














## Grazing model
BehSum5$Not_graze <- BehSum5$Total - BehSum5$Graze
Graze_mod <- glmmTMB(cbind(Graze, Not_graze) ~ pre_post + (1|Tag_ID:pair_ID),
                     data = BehSum5,
                     family = betabinomial)
summary(Graze_mod)
drop1(Graze_mod, test = "Chi")
model_performance(Graze_mod) 
check_model(Graze_mod) 


## Stationary model 
BehSum5$Not_Stat <- BehSum5$Total - BehSum5$Stat
Stat_mod <- glmmTMB(cbind(Stat, Not_Stat) ~ pre_post + (1|Tag_ID:pair_ID),
                     data = BehSum5,
                     family = betabinomial)
summary(Stat_mod)
drop1(Stat_mod, test = "Chi")
model_performance(Stat_mod) 
check_model(Stat_mod) 


## Flight model 
BehSum5$Not_Flight <- BehSum5$Total - BehSum5$Flight
Flight_mod <- glmmTMB(cbind(Flight, Not_Flight) ~ pre_post + (1|Tag_ID:pair_ID),
                    data = BehSum5,
                    family = betabinomial)
summary(Flight_mod)
drop1(Flight_mod, test = "Chi")
model_performance(Flight_mod) 
check_model(Flight_mod) 


## Flight model 
BehSum5$Not_Walk <- BehSum5$Total - BehSum5$Walk
Walk_mod <- glmmTMB(cbind(Walk, Not_Walk) ~ pre_post + (1|Tag_ID:pair_ID),
                      data = BehSum5,
                      family = betabinomial)
summary(Walk_mod)
drop1(Walk_mod, test = "Chi")
model_performance(Walk_mod) 
check_model(Walk_mod) 



