## Luke Ozsanlav-Harris
## Created 29/04/2022

## Collate together ODBA data from GBG and GWfG
## Test the effects on different burst lengths on daily ODBA 
## then model if ODBA changes on shooting and non-shooting days


## ## Packages required
pacman::p_load(tidyverse, lubridate, data.table, zoo, svMisc, suncalc,
               geosphere, sf, sfheaders, DHARMa, MuMIn, effects, 
               performance, glmmTMB, emmeans, raster)



#------------------------------------------#
#### 1. Read in the GWfG ODBA data sets ####
#------------------------------------------#

#-----------------------------------#
#### 1.1 GWfG Ornitela ODBA data ####
#-----------------------------------#

## read in the quantile mapped ODBA from the ornitela tags
GWfG_Orn <- readRDS("tracking data/All_orni_ODBAQuantile.RDS")


## parse timestamp
GWfG_Orn$UTC_datetime <- ymd_hms(GWfG_Orn$UTC_datetime)


## Check for duplicated observations (ones with same lat, long, timestamp, and device ID)
ind2 <- GWfG_Orn %>% 
        dplyr::select(UTC_datetime, device_id) %>%
        duplicated()

GWfG_Orn <- GWfG_Orn %>% filter(!ind2)

## Now filter for the winter, create a year day column first
GWfG_Orn <- GWfG_Orn %>% 
            mutate(year_day = yday(UTC_datetime)) %>% 
            filter(year_day >= 305 | year_day <= 92)




#----------------------------------#
#### 1.2 GWfG Ecotone ODBA data ####
#----------------------------------#

## read in the Ecotone ODBA
GWfG_Eco <- readRDS("tracking data/Eco_ODBA_per_burst.RDS")

## parse timestamp
GWfG_Eco$UTC_datetime <- ymd_hms(GWfG_Eco$UTC_datetime)

##Check for duplicated observations (ones with same lat, long, timestamp, and device ID)
ind3 <- GWfG_Eco %>% 
        dplyr::select(UTC_datetime, Tag_ID) %>%
        duplicated()

GWfG_Eco <- GWfG_Eco %>% filter(!ind3)

## Now filter for the winter, create a year day column first
GWfG_Eco <- GWfG_Eco %>% 
            mutate(year_day = yday(UTC_datetime)) %>% 
            filter(year_day >= 305 | year_day <= 92)




#-------------------------------------------#
#### 1.3 Bind two GWfG datasets together ####
#-------------------------------------------#

## Streamline the two datasets
GWfG_Orn <- subset(GWfG_Orn, select= c("device_id", "UTC_datetime", "ODBA", "ODBA_quant", "year_day"))
GWfG_Eco <- subset(GWfG_Eco, select= c("Tag_ID", "UTC_datetime", "ODBA", "year_day"))
setnames(GWfG_Eco, old = "Tag_ID", new = "device_id")

## Bind the two datasets together for GWfG
winter_GWfG <- plyr::rbind.fill(GWfG_Orn, GWfG_Eco)


## Add sex data to the GWfG data
Add_data <- fread("MetaData/Tagged bird summary data new.csv")

## sort out bird IDs
Add_data <- Add_data %>% 
            mutate(S.N= as.character(Add_data$S.N),
                   S.N= ifelse(is.na(S.N)==T, Bird.ID, S.N)) %>% 
            dplyr::select(S.N, Sex)

## Join sexes to data set
ind4 <- duplicated(Add_data)
Add_data <- filter(Add_data, !ind4)
winter_GWfG2 <- left_join(winter_GWfG, Add_data, by =c("device_id"="S.N"))




#---------------------------------------------#
#### 2. Read in GWfG Islay winter GPS data ####
#---------------------------------------------#

## read in GPS data
GPS_data <- readRDS("Derived data/All_winter_GPS_with_habitat.RDS")

## give each fix its own unique ID number
GPS_data$ID <- 1:nrow(GPS_data)




#-----------------------------------------------------------#
#### 3. Subset GWfG ACC data to just the extent of Islay ####
#-----------------------------------------------------------#

## The GPS_set data set contains just GPS data from the winter on Islay
## Therefore if i restrict the ACC data to the same tag periods then i can constrain the ACC data to just come from Islay winter

## Add a tag date column to both data sets
GPS_data$tag_date <- paste0(GPS_data$Tag_ID, "_", as.Date(GPS_data$UTC_datetime))
winter_GWfG2$tag_date <- paste0(winter_GWfG2$device_id, "_", as.Date(winter_GWfG2$UTC_datetime))

## Add a tag winter column to both data sets
GPS_data$winter <- ifelse(GPS_data$year_day < 150, paste0((year(GPS_data$UTC_datetime)-1), "-", year(GPS_data$UTC_datetime)), paste0(year(GPS_data$UTC_datetime), "-", (year(GPS_data$UTC_datetime)+1)))
winter_GWfG2$winter <- ifelse(winter_GWfG2$year_day < 150, paste0((year(winter_GWfG2$UTC_datetime)-1), "-", year(winter_GWfG2$UTC_datetime)), paste0(year(winter_GWfG2$UTC_datetime), "-", (year(winter_GWfG2$UTC_datetime)+1)))

## Now filter out only the tag dates in the GPS data set
UnQs <- unique(GPS_data$tag_date)
winter_GWfG2 <- filter(winter_GWfG2, tag_date %in% UnQs)


## Now filter out bird winter when there was 5 or less days of ODBA data 
## first create summary of the number of days of ODBA data for each tag winter
ODBA_sum <- winter_GWfG2 %>% 
            group_by(device_id, winter) %>% 
            summarise(n_days = length(unique(tag_date)))

## now filter out instances with 5 or less days
ODBA_min <- filter(ODBA_sum, n_days > 5)
ODBA_min$n_days <- NULL

## Use inner join to remove the tag winters with 5 or less days
winter_GWfG3 <- inner_join(winter_GWfG2, ODBA_min, by = c("device_id", "winter"))

## remove tag 17795 as the ACC sensor appear to be non-functional at times
winter_GWfG3 <- filter(winter_GWfG3, !device_id == 17795)




#------------------------------------#
#### 4. Read in the GBG ODBA data ####
#------------------------------------#

## read in ODBA data
GBG <- fread("tracking data/ALL_TAGS_ODBA_GPS.csv")
GBG$V1 <- NULL

## parse timestamp
GBG$Timestamp_ODBA <- ymd_hms(GBG$Timestamp_ODBA)

##Check for duplicated observations (ones with same lat, long, timestamp, and device ID)
ind <- GBG %>% 
       dplyr::select(Timestamp_ODBA, device_id) %>%
       duplicated()

GBG <- GBG %>% filter(!ind)

## Now filter for the winter, create a year day column first
GBG <- GBG %>% 
       mutate(year_day = yday(Timestamp_ODBA)) %>% 
       filter(year_day >= 305 | year_day <= 92)


## trip the data to just the extent of Islay
## Read in a spatial object with the extent of Islay
LandCov_2015 <- raster("Landcover Data/Islay landcover data/LandCov2015_Islay.grd")

## Convert the GPS data into a sf object
GBG_sf <- st_as_sf(GBG, coords = c("Longitude", "Latitude"), 
                    crs = 4326, agr = "constant")

## now convert to the same crs as the field data
GBG_sf <- st_transform(GBG_sf, crs = st_crs(LandCov_2015))

## crop to the same extent as the landcover data
Islay_GBG_sf <- st_crop(GBG_sf, st_bbox(LandCov_2015))

## covert the sf object back to a data frame
GBG_ODBA <- sf_to_df(Islay_GBG_sf, fill = T)




#-----------------------------------------------#
#### 5. Join together GBG and GWfG ODBA data ####
#-----------------------------------------------#

## subset the GWfG data
winter_GWfG3 <- winter_GWfG3 %>% 
                mutate(ODBA_quant = ifelse(is.na(ODBA_quant) == T, ODBA, ODBA_quant)) %>% 
                dplyr::select(device_id, UTC_datetime, Sex, ODBA_quant, year_day) %>% 
                rename(ODBA= ODBA_quant)

## subset GBG data set
GBG_ODBA <- subset(GBG_ODBA, select = c("device_id", "Timestamp_ODBA", "Sex", "ODBA", "year_day"))
setnames(GBG_ODBA, old = "Timestamp_ODBA", new = "UTC_datetime")

## fix sex column in GBG
GBG_ODBA$Sex <- ifelse(GBG_ODBA$Sex == "MALE", "M", "F")

## set species column
GBG_ODBA$Species <- "GBG"
winter_GWfG3$Species <- "GWfG"

## bind the two together
All_winter_ODBA <- plyr::rbind.fill(winter_GWfG3, GBG_ODBA)
table(All_winter_ODBA$device_id); unique(All_winter_ODBA$Sex)




#------------------------------------------#
#### 6. Identify days birds are Shot on ####
#------------------------------------------#

#-----------------------------------------------------------#                                               
#### 6.1 Read in GPS tracking and join two species data ####
#-----------------------------------------------------------#  

## read in GWfG GPS data 
GWfG_GPS <- readRDS("Derived data/All_winter_GPS_with_habitat.RDS")

## read in the tracking data csv for GBG
GBG_GPS <- fread("tracking data/GBG_TAG_CLEAN.csv")

## extract columns that i want from the GBG data
GBG_GPS_sub <- subset(GBG_GPS, select = c("device_id", "UTC_datetime", "Longitude", "Latitude"))

## change column names of the GBG data so it matches the GWfG data, add year day column and parse timestamp first
GBG_GPS_sub <- GBG_GPS_sub %>% 
               rename(Tag_ID = device_id) %>% 
               mutate(species = "GBG",
                      UTC_datetime = dmy_hms(UTC_datetime), 
                      year_day = yday(UTC_datetime)) %>% 
               filter(year_day >= 305 | year_day <= 92)

## Now join to the GWfG data, make sure the timestamps are the same class
GWfG_GPS$UTC_datetime <- ymd_hms(GWfG_GPS$UTC_datetime)
GBG_GPS_sub$UTC_datetime <- ymd_hms(GBG_GPS_sub$UTC_datetime)
GWfG_GPS$Habitat <- NULL
winter_gps  <- rbind(GWfG_GPS, GBG_GPS_sub) 

##  order by TagID and by timestamp then add a unique ID to each fix
winter_gps <- winter_gps[order(winter_gps$Tag_ID, winter_gps$UTC_datetime),]
winter_gps$ID <- 1:nrow(winter_gps)




#-------------------------------------------------------#                                                
#### 6.2 Read in shooting data and field boundaries ####
#-------------------------------------------------------#                                               

## read in the shooting logs
Logs <- fread("Shooting logs/All_logs_cleaned.csv")

## create a timestamp column for the shooting logs, think the ones that fail to parse are the ones were the time is missing and it is just a date
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




#--------------------------------------------------------#                                               
#### 6.3 Combine shooting data and spatial field data ####
#--------------------------------------------------------#                                                    

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




#-----------------------------------------------------------#                                                   
#### 6.4 Label fixes associated with each shooting event ####
#-----------------------------------------------------------#                                                     

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
  for (i in 1:nrid) {
    
    ## progress indicator for the loop
    svMisc::progress(value = i, max.value = nrid)
    
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


## run the function on the species separately
## first split the species
GBG_sub <- filter(winter_gps, species == "GBG")
GWfG_sub <- filter(winter_gps, species == "GWfG")

## run the function on both data sets with different thresholds
Prox_list1 <- Shoot_proximity(TD = GBG_sub, SD = winter_shots_cent2, time_tresh = 60, dist_thresh = 1184, field_crs = Field_centres)
Prox_list2 <- Shoot_proximity(TD = GWfG_sub, SD = winter_shots_cent2, time_tresh = 60, dist_thresh = 644, field_crs = Field_centres)

## bind the lists together into one data frame
All_prox1 <- plyr::ldply(Prox_list1)
All_prox2 <- plyr::ldply(Prox_list2)

## Join the species lists
All_prox <- rbind(All_prox1, All_prox2)


#-----------------------------------------------------------------------#                                                                    
#### 7. Summarize on what days and how many times birds were shot at ####
#-----------------------------------------------------------------------#                                                                  

## create a date column for the grouping
All_prox$date <- as.Date(All_prox$UTC_datetime)

## create a dummy column to sum rows
All_prox$dummy <- 1

## now summarise the  number of shooting events by each tag day
DistSum <- All_prox %>% 
            group_by(date, Tag_ID) %>% 
            summarise(n_shot = sum(dummy))




#---------------------------------------#                                                                   
#### 8. Summarize average daily ODBA ####
#---------------------------------------# 

## First label times as either day or night
## "sunrise" and "sunset" gave short days, i.e. the geese would still be out on daytime feeding sites
## Either use  "dusk", "dawn" or "nauticalDusk", "nauticalDawn" (nautical versions give the longer day length)
sun_dates <- seq.Date(as.Date(min(All_winter_ODBA$UTC_datetime)), as.Date(max(All_winter_ODBA$UTC_datetime)), by = 1)
sunt <-
  getSunlightTimes(
    date = sun_dates,
    keep = c("dusk", "dawn"),
    lat = 55.738,
    lon = 6.246,
    tz = "UTC"
  )


## Add a column for the date
All_winter_ODBA$date <- as.Date(All_winter_ODBA$UTC_datetime)

## Now join on the susnet and sunrise times
sunt$lat <- NULL
sunt$lon <- NULL
All_winter_ODBA2 <- left_join(All_winter_ODBA, sunt,  by = "date")
All_winter_ODBA2$DayNight <- ifelse(All_winter_ODBA2$UTC_datetime < All_winter_ODBA2$dawn | All_winter_ODBA2$UTC_datetime > All_winter_ODBA2$dusk,
                                    "night", "day")


## filter out ODBA values of 0
All_winter_ODBA2 <- filter(All_winter_ODBA2, !ODBA == 0)


## Create summary of behaviors before and after shooting
All_winter_ODBA2$dummy <- 1
ODBASum <- All_winter_ODBA2 %>% 
           group_by(date, device_id, DayNight, Sex, Species) %>% 
           summarise(Avg_ODBA = mean(ODBA),
                     Var_ODBA = var(ODBA),
                     Sum_ODBA = sum(ODBA),
                     n_readings = sum(dummy),
                     StEr_ODBA = sqrt(Var_ODBA / n_readings))




#-------------------------------------------------------# 
#### 9. Add on number of shooting events to the data ####
#-------------------------------------------------------# 

## join number of shots and daily behavioral summary
setnames(ODBASum, old = "device_id", new = "Tag_ID")
ODBASum2 <- left_join(ODBASum, DistSum, by = c("date", "Tag_ID"))
stopifnot(nrow(ODBASum2)==nrow(ODBASum2))


## change all the NAs in the n_shot column to 0s
ODBASum2$n_shot <- ifelse(is.na(ODBASum2$n_shot)==T, 0, ODBASum2$n_shot)




#---------------------------------------------------------------------# 
#### 10. Model Average daily ODBA on shooting vs non-shooting days ####
#---------------------------------------------------------------------#

#-----------------------------------------#
#### 10.1 Add columns for use in model ####
#-----------------------------------------#

## Add extra columns to put in to the model
ODBASum2$year <- year(ODBASum2$date) # year column
ODBASum2$year_day <- yday(ODBASum2$date) 
ODBASum2$from_nov1 <- ifelse(ODBASum2$year_day > 250, 
                            ODBASum2$year_day - 305, 
                            ODBASum2$year_day + 60) # day since November 1st
ODBASum2$from_nov1Factor <- as.factor(ODBASum2$from_nov1)
ODBASum2$shot <- ifelse(ODBASum2$n_shot >  0, "Shot", "Not_Shot") # whether there was shooting or not
ODBASum2$winter <- ifelse(ODBASum2$year_day < 150, paste0((ODBASum2$year-1), "-", ODBASum2$year), paste0(ODBASum2$year, "-", (ODBASum2$year+1)))
ODBASum2$tag_winter <- paste0(ODBASum2$Tag_ID, "_", ODBASum2$winter) # tag winter identifier
ODBASum2$prev_shot <- ifelse(lag(ODBASum2$shot, n = 2L) == "Shot", "Prev_shot", "Not")
ODBASum2$prev_shot <- ifelse(is.na(ODBASum2$prev_shot)==T, "Not", ODBASum2$prev_shot)

## set variables as correct class
ODBASum2 <- ODBASum2 %>% ungroup() %>%  mutate(across(c(year, DayNight, Tag_ID, shot, from_nov1), as.factor))
ODBASum2$from_nov1 <- as.numeric(ODBASum2$from_nov1)

## Create another column for the inverse of the varience and standard errors
ODBASum2$InVar_ODBA <- 1/ODBASum2$Var_ODBA
ODBASum2$InStEr_ODBA <- 1/ODBASum2$StEr_ODBA




#-------------------------#
#### 10.2 Prelim plots ####
#-------------------------#

## plot by species and Tag ID
ggplot(data = ODBASum2) +
  geom_boxplot(aes(Avg_ODBA, Tag_ID, colour = Species)) + 
  theme_bw()

## plot by GWfG Tag ID
ODBASum2 %>% 
  filter(Species == "GWfG") %>% 
  ggplot() +
    geom_boxplot(aes(Avg_ODBA, Tag_ID)) + 
    theme_bw()

## plot by GWfG Tag ID
ODBASum2 %>% 
  filter(Species == "GBG") %>% 
  ggplot() +
  geom_boxplot(aes(Avg_ODBA, Tag_ID)) + 
  theme_bw()

## historgam of Avg_ODBA
hist(sqrt(ODBASum2$Avg_ODBA))




#--------------------------------------------#
#### 10.3 Minimum data requirement checks ####
#--------------------------------------------#

## remove data points based on whether they had a certain number of data points per day
ODBASum3 <- filter(ODBASum2, n_readings >= 6)

## now filter out instances where an individual has 5 or less day/night worth of data
ID_sum <- table(ODBASum3$tag_winter) %>% data.frame()
ID_sum <- filter(ID_sum, Freq > 10)


## Use inner join to remove the tag winters with 5 or less days
ODBASum4 <- filter(ODBASum3, tag_winter %in% ID_sum$Var1)

## create summary to see how much data I have
datasize <- ODBASum4 %>% group_by(Species, DayNight, n_shot) %>% summarise(samp = n())
datasize

## Low samples size of shooting disturbed at night for GBG, so filtering out night for both species
ODBASum4 <- filter(ODBASum4, DayNight == "day")




#--------------------------------------------------------------#
#### 10.4 Test multiple models with different terms in them ####
#--------------------------------------------------------------#

## Start with the interaction term but then gradually drop the terms and see which has the lowest AIC

Mod_Int <-  glmmTMB(formula = Sum_ODBA ~ n_readings + shot*Species + Sex + (1|winter) + poly(from_nov1,2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                   data = ODBASum4,
                   family = Gamma(link = "log"),
                   REML = FALSE,
                   control=glmmTMBControl(optimizer=optim,
                                          optArgs=list(method="BFGS")))

Mod_NoInt <-  glmmTMB(formula = Sum_ODBA ~ n_readings + shot + Species + Sex + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                   data = ODBASum4,
                   family = Gamma(link = "log"),
                   REML = FALSE,
                   control=glmmTMBControl(optimizer=optim,
                                          optArgs=list(method="BFGS")))

Mod_Species <-  glmmTMB(formula = Sum_ODBA ~ n_readings + Species + Sex + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                   data = ODBASum4,
                   family = Gamma(link = "log"),
                   REML = FALSE,
                   control=glmmTMBControl(optimizer=optim,
                                          optArgs=list(method="BFGS")))

Mod_Shot <-  glmmTMB(formula = Sum_ODBA ~ n_readings + shot + Sex + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                      data = ODBASum4,
                      family = Gamma(link = "log"),
  
                      REML = FALSE,
                      control=glmmTMBControl(optimizer=optim,
                                             optArgs=list(method="BFGS")))

Mod_None <-  glmmTMB(formula = Sum_ODBA ~ n_readings + Sex + (1|winter) + poly(from_nov1, 2) + (1|Tag_ID) + ar1(from_nov1Factor + 0 | tag_winter),
                   data = ODBASum4,
                   family = Gamma(link = "log"),
                   REML = FALSE,
                   control=glmmTMBControl(optimizer=optim,
                                          optArgs=list(method="BFGS")))

## compare models with AIC
AIC(Mod_Int, Mod_NoInt, Mod_Species, Mod_Shot, Mod_None)

## Now get the estimates from the main mode
summary(Mod_None) # quick summary
confint(Mod_None)[1:10,] #confidence intervals
MuMIn::r.squaredGLMM(Mod_None)
confint(Mod_Int)

## model checks
ResidsMod_ReadingsWeight <- simulateResiduals(Mod_Int)
plot(ResidsMod_ReadingsWeight)
par(mfrow=c(1,3))
plotResiduals(ResidsMod_ReadingsWeight[["scaledResiduals"]], form = ODBASum4$DayNight)
plotResiduals(ResidsMod_ReadingsWeight[["scaledResiduals"]], form = ODBASum4$shot)
plotResiduals(ResidsMod_ReadingsWeight[["scaledResiduals"]], form = ODBASum4$from_nov1)


## the below computes all of the pairwise comparisons
## and then calculates the confidence intervals for these
## not sure this is the right way of doing this.....
## create all the pairwise comparisons within the two way interaction
grid1 <- emmeans(Mod_Int, ~shot*Species)
PairComp <- as.data.frame(pairs(grid1))

## Now estimate the CIs for each of the estimates
PairComp$lowCI <- PairComp$estimate - (1.96*PairComp$SE)
PairComp$highCI <- PairComp$estimate + (1.96*PairComp$SE)




#------------------------------------------------#
#### 10.5 Plot the output from the best model ####
#------------------------------------------------#

## use the effects package, and ggplot to plot the model predictions
## first the effects for each predicitor
top_mod_effects <- predictorEffects(Mod_None)
plot(top_mod_effects)

## now extract the fits for the first variable and bind them together
effects3 <- top_mod_effects["shot"]
fit3 <- as.data.frame(cbind(as.numeric(effects3[["shot"]][["fit"]]), as.numeric(effects3[["shot"]][["lower"]]), 
                            as.numeric(effects3[["shot"]][["upper"]]), as.character(effects3[["shot"]][["x"]][["shot"]]),
                            as.character(effects3[["shot"]][["x"]][["Species"]])))
## change the names to something meaningful
setnames(fit3, old = c("V1", "V2", "V3", "V4", "V5"), new = c("fit", "lower", "upper", "shot", "Species"))
## transform variables back to proportion scale
fit3$fit <- exp(as.numeric(as.character(fit3$fit)))
fit3$lower <- exp(as.numeric(as.character(fit3$lower)))
fit3$upper <- exp(as.numeric(as.character(fit3$upper)))

## Change the names of some of the variables for plotting
fit3$shot <- ifelse(fit3$shot == "Shot", "Shooting", "No Shooting")
ODBASum5 <- ODBASum4
ODBASum5$shot <- ifelse(ODBASum5$shot == "Shot", "Shooting", "No Shooting")

## Now plot using ggplot
ggplot() + 
  geom_errorbar(data = fit3, aes(x= shot, ymin = lower, ymax = upper, colour = Species), width = 0.4, size = 0.8) +
  geom_point(data=fit3, aes(x= shot, y = fit, colour = Species), size = 2, shape = 15)  +
  #geom_point(data = ODBASum5, aes(x =shot, y= Sum_ODBA, colour = Species), size = 0.75, alpha = 0.1) +
  xlab("Shooting Exposure") + ylab("Average daily ODBA") +
  theme_bw() +
  facet_wrap(~Species)+
  scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  labs(colour="Species") +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=12,), 
        panel.grid.minor.x = element_blank(), 
        legend.position = "none",
        panel.grid.major.x = element_blank(), 
        axis.text = element_text(size =14), 
        axis.title.x = element_text(size =18),
        axis.title.y = element_text(size =18), 
        strip.text.x = element_text(size =14), 
        legend.title = element_text(size =14),
        legend.text = element_text(size =14))

## Save a plot
# ggsave("Plots/Script 4) plots/Average daily ODBA by species.png", 
#        width = 18, height = 20, units = "cm")

